#include <algorithm>

#include "NetworkInfo.hpp"
#include "ParallelContext.hpp"

using namespace std;

void NetworkInfo::init(const Options &opts, const Network& network, const PartitionedMSA& parted_msa,
                    const IDVector& tip_msa_idmap,
                    const PartitionAssignment& part_assign,
                    const std::vector<uintVector>& site_weights)
{
  _brlen_min = opts.brlen_min;
  _brlen_max = opts.brlen_max;
  _brlen_opt_method = opts.brlen_opt_method;
  _partition_contributions.resize(parted_msa.part_count());
  double total_weight = 0;

  _pll_networkinfo = pllmod_networkinfo_create(pll_unetwork_clone(&network.pll_unetwork()),
                                         network.num_tips(),
                                         parted_msa.part_count(), opts.brlen_linkage);

  libpll_check_error("ERROR creating networkinfo structure");
  assert(_pll_networkinfo);

  if (ParallelContext::num_procs() > 1)
  {
    pllmod_networkinfo_set_parallel_context(_pll_networkinfo, (void *) nullptr,
                                         ParallelContext::parallel_reduce_cb);
  }

  // init partitions
  int optimize_branches = opts.optimize_brlen ? PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE : 0;

  for (size_t p = 0; p < parted_msa.part_count(); ++p)
  {
    const PartitionInfo& pinfo = parted_msa.part_info(p);
    const auto& weights = site_weights.empty() ? pinfo.msa().weights() : site_weights.at(p);
    int params_to_optimize = opts.optimize_model ? pinfo.model().params_to_optimize() : 0;
    params_to_optimize |= optimize_branches;

    _partition_contributions[p] = std::accumulate(weights.begin(), weights.end(), 0);
    total_weight += _partition_contributions[p];

    PartitionAssignment::const_iterator part_range = part_assign.find(p);
    if (part_range != part_assign.end())
    {
      /* create and init PLL partition structure */
      pll_partition_t * partition = create_pll_partition_network(opts, pinfo, tip_msa_idmap,
                                                         *part_range, weights);

      int retval = pllmod_networkinfo_init_partition(_pll_networkinfo, p, partition,
                                                  params_to_optimize,
                                                  pinfo.model().gamma_mode(),
                                                  pinfo.model().alpha(),
                                                  pinfo.model().ratecat_submodels().data(),
                                                  pinfo.model().submodel(0).rate_sym().data());

      if (!retval)
      {
        assert(pll_errno);
        libpll_check_error("ERROR adding treeinfo partition");
      }

      // set per-partition branch lengths or scalers
      if (opts.brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
      {
        assert (_pll_networkinfo->brlen_scalers);
        _pll_networkinfo->brlen_scalers[p] = pinfo.model().brlen_scaler();
      }
      else if (opts.brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED && !network.partition_brlens().empty())
      {
        assert(_pll_networkinfo->branch_lengths[p]);
        memcpy(_pll_networkinfo->branch_lengths[p], network.partition_brlens(p).data(),
               network.num_branches() * sizeof(double));
      }

      if (part_range->master())
        _parts_master.insert(p);
    }
    else
    {
      // this partition will be processed by other threads, but we still need to know
      // which parameters to optimize
      _pll_networkinfo->params_to_optimize[p] = params_to_optimize;
    }
  }

  // finalize partition contribution computation
  for (auto& c: _partition_contributions)
    c /= total_weight;
}

NetworkInfo::NetworkInfo (const Options &opts, const Network& network, const PartitionedMSA& parted_msa,
                    const IDVector& tip_msa_idmap,
                    const PartitionAssignment& part_assign) : _pll_networkinfo(nullptr), _brlen_max(0), _brlen_min(0), _brlen_opt_method(0)
{
  init(opts, network, parted_msa, tip_msa_idmap, part_assign, std::vector<uintVector>());
}

NetworkInfo::NetworkInfo (const Options &opts, const Network& network, const PartitionedMSA& parted_msa,
                    const IDVector& tip_msa_idmap,
                    const PartitionAssignment& part_assign,
                    const std::vector<uintVector>& site_weights) : _pll_networkinfo(nullptr), _brlen_max(0), _brlen_min(0), _brlen_opt_method(0)
{
  init(opts, network, parted_msa, tip_msa_idmap, part_assign, site_weights);
}

NetworkInfo::~NetworkInfo ()
{
  if (_pll_networkinfo)
  {
    for (unsigned int i = 0; i < _pll_networkinfo->partition_count; ++i)
    {
      if (_pll_networkinfo->partitions[i])
        pll_partition_destroy(_pll_networkinfo->partitions[i]);
    }

    pll_unetwork_graph_destroy(_pll_networkinfo->root, NULL);
    pllmod_networkinfo_destroy(_pll_networkinfo);
  }
}

Network NetworkInfo::network() const
{
  if (!_pll_networkinfo)
    return Network();

  Network tree(*_pll_networkinfo->network);

  if (_pll_networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    // set per-partition branch lengths
    for (unsigned int i = 0; i < _pll_networkinfo->partition_count; ++i)
    {
      assert(_pll_networkinfo->branch_lengths[i]);
      doubleVector brlens(_pll_networkinfo->branch_lengths[i],
    		  _pll_networkinfo->branch_lengths[i] + tree.num_branches());
      tree.add_partition_brlens(std::move(brlens));
    }

    // compute a weighted average of per-partition brlens
    tree.apply_avg_brlens(_partition_contributions);
  }

  return tree;
}

Network NetworkInfo::network(size_t partition_id) const
{
  if (!_pll_networkinfo)
    return Network();

  if (partition_id >= _pll_networkinfo->partition_count)
    throw out_of_range("Partition ID out of range");

  PllUNetworkUniquePtr pll_unetwork(pllmod_networkinfo_get_partition_network(_pll_networkinfo, partition_id));

  if (!pll_unetwork)
  {
    assert(pll_errno);
    libpll_check_error("networkinfo: cannot get partition network");
  }

  return Network(pll_unetwork);
}

void NetworkInfo::network(const Network& network)
{
	_pll_networkinfo->network = pll_unetwork_clone(&network.pll_unetwork());
	_pll_networkinfo->root = _pll_networkinfo->network->vroot;
}

double NetworkInfo::loglh(bool incremental)
{
  return pllmod_networkinfo_compute_loglh(_pll_networkinfo, incremental ? 1 : 0, 1);
}

void NetworkInfo::model(size_t partition_id, const Model& model)
{
  if (partition_id >= _pll_networkinfo->partition_count)
    throw out_of_range("Partition ID out of range");

  if (!_pll_networkinfo->partitions[partition_id])
    return;

  assign(_pll_networkinfo->partitions[partition_id], model);
  _pll_networkinfo->alphas[partition_id] = model.alpha();
  if (_pll_networkinfo->brlen_scalers)
	  _pll_networkinfo->brlen_scalers[partition_id] = model.brlen_scaler();
}

//#define DBG printf

double NetworkInfo::optimize_branches(double lh_epsilon, double brlen_smooth_factor)
{
  /* update all CLVs and p-matrices before calling BLO */
  double new_loglh = loglh();

  if (_pll_networkinfo->params_to_optimize[0] & PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE)
  {
    int max_iters = brlen_smooth_factor * RAXML_BRLEN_SMOOTHINGS;
    new_loglh = -1 * pllmod_algo_opt_brlen_networkinfo(_pll_networkinfo,
                                                    _brlen_min,
                                                    _brlen_max,
                                                    lh_epsilon,
                                                    max_iters,
                                                    _brlen_opt_method,
                                                    PLLMOD_OPT_BRLEN_OPTIMIZE_ALL
                                                    );

    LOG_DEBUG << "\t - after brlen: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in branch length optimization");
    assert(isfinite(new_loglh));
  }

  /* optimize brlen scalers, if needed */
  if (_pll_networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED &&
		  _pll_networkinfo->partition_count > 1)
  {
    new_loglh = -1 * pllmod_algo_opt_brlen_scalers_networkinfo(_pll_networkinfo,
                                                            RAXML_BRLEN_SCALER_MIN,
                                                            RAXML_BRLEN_SCALER_MAX,
                                                            _brlen_min,
                                                            _brlen_max,
                                                            RAXML_PARAM_EPSILON);

    LOG_DEBUG << "\t - after brlen scalers: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in brlen scaler optimization");
    assert(isfinite(new_loglh));
  }

  return new_loglh;
}

double NetworkInfo::optimize_params(int params_to_optimize, double lh_epsilon)
{
  assert(!pll_errno);

  double
    cur_loglh = loglh(),
    new_loglh = cur_loglh;

  /* optimize SUBSTITUTION RATES */
  if (params_to_optimize & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_subst_rates_networkinfo(_pll_networkinfo,
                                                          0,
                                                          PLLMOD_OPT_MIN_SUBST_RATE,
                                                          PLLMOD_OPT_MAX_SUBST_RATE,
                                                          RAXML_BFGS_FACTOR,
                                                          RAXML_PARAM_EPSILON);

    LOG_DEBUG << "\t - after rates: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in substitution rates optimization");
    assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
    cur_loglh = new_loglh;
  }

  /* optimize BASE FREQS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREQUENCIES)
  {
    new_loglh = -1 * pllmod_algo_opt_frequencies_networkinfo(_pll_networkinfo,
                                                          0,
                                                          PLLMOD_OPT_MIN_FREQ,
                                                          PLLMOD_OPT_MAX_FREQ,
                                                          RAXML_BFGS_FACTOR,
                                                          RAXML_PARAM_EPSILON);

    LOG_DEBUG << "\t - after freqs: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in base frequencies optimization");
    assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
    cur_loglh = new_loglh;
  }

  // TODO: co-optimization of PINV and ALPHA, mb with multiple starting points
  if (0 &&
      (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA) &&
      (params_to_optimize & PLLMOD_OPT_PARAM_PINV))
  {
    new_loglh = -1 * pllmod_algo_opt_alpha_pinv_networkinfo(_pll_networkinfo,
                                                         0,
                                                         PLLMOD_OPT_MIN_ALPHA,
                                                         PLLMOD_OPT_MAX_ALPHA,
                                                         PLLMOD_OPT_MIN_PINV,
                                                         PLLMOD_OPT_MAX_PINV,
                                                         RAXML_BFGS_FACTOR,
                                                         RAXML_PARAM_EPSILON);

    LOG_DEBUG << "\t - after a+i  : logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in alpha/p-inv parameter optimization");
    assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
    cur_loglh = new_loglh;
  }
  else
  {
    /* optimize ALPHA */
    if (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA)
    {
      new_loglh = -1 * pllmod_algo_opt_onedim_networkinfo(_pll_networkinfo,
                                                        PLLMOD_OPT_PARAM_ALPHA,
                                                        PLLMOD_OPT_MIN_ALPHA,
                                                        PLLMOD_OPT_MAX_ALPHA,
                                                        RAXML_PARAM_EPSILON);

     LOG_DEBUG << "\t - after alpha: logLH = " << new_loglh << endl;

     libpll_check_error("ERROR in alpha parameter optimization");
     assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
     cur_loglh = new_loglh;
    }

    /* optimize PINV */
    if (params_to_optimize & PLLMOD_OPT_PARAM_PINV)
    {
      new_loglh = -1 * pllmod_algo_opt_onedim_networkinfo(_pll_networkinfo,
                                                        PLLMOD_OPT_PARAM_PINV,
                                                        PLLMOD_OPT_MIN_PINV,
                                                        PLLMOD_OPT_MAX_PINV,
                                                        RAXML_PARAM_EPSILON);

      LOG_DEBUG << "\t - after p-inv: logLH = " << new_loglh << endl;

      libpll_check_error("ERROR in p-inv optimization");
      assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
      cur_loglh = new_loglh;
    }
  }

  /* optimize FREE RATES and WEIGHTS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREE_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_rates_weights_networkinfo (_pll_networkinfo,
                                                          RAXML_FREERATE_MIN,
                                                          RAXML_FREERATE_MAX,
                                                          RAXML_BFGS_FACTOR,
                                                          RAXML_PARAM_EPSILON);

    /* normalize scalers and scale the branches accordingly */
    if (_pll_networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED &&
    		_pll_networkinfo->partition_count > 1)
      pllmod_networkinfo_normalize_brlen_scalers(_pll_networkinfo);

    LOG_DEBUG << "\t - after freeR: logLH = " << new_loglh << endl;
//    LOG_DEBUG << "\t - after freeR/crosscheck: logLH = " << loglh() << endl;

    libpll_check_error("ERROR in FreeRate rates/weights optimization");
    assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
    cur_loglh = new_loglh;
  }

  if (params_to_optimize & PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE)
  {
    new_loglh = optimize_branches(lh_epsilon, 0.25);

    assert(cur_loglh - new_loglh < -new_loglh * RAXML_DOUBLE_TOLERANCE);
    cur_loglh = new_loglh;
  }

  return new_loglh;
}

double NetworkInfo::spr_round(spr_round_params& params)
{
  double loglh = pllmod_algo_spr_round_network(_pll_networkinfo, params.radius_min, params.radius_max,
                               params.ntopol_keep, params.thorough, _brlen_opt_method,
                               _brlen_min, _brlen_max, RAXML_BRLEN_SMOOTHINGS,
                               0.1,
                               params.subtree_cutoff > 0. ? &params.cutoff_info : nullptr,
                               params.subtree_cutoff);

  libpll_check_error("ERROR in SPR round");

  assert(isfinite(loglh) && loglh);

  return loglh;
}

void NetworkInfo::set_topology_constraint(const Network& cons_network)
{
  if (!cons_network.empty())
  {
    int retval = pllmod_networkinfo_set_constraint_network(_pll_networkinfo, &cons_network.pll_unetwork());
    if (!retval)
      libpll_check_error("ERROR: Cannot set topological constraint");
  }
}

void assign(PartitionedMSA& parted_msa, const NetworkInfo& networkinfo)
{
  const pllmod_networkinfo_t& pll_networkinfo = networkinfo.pll_networkinfo();

  if (parted_msa.part_count() != pll_networkinfo.partition_count)
    throw runtime_error("Incompatible arguments");

  for (size_t p = 0; p < parted_msa.part_count(); ++p)
  {
    if (!pll_networkinfo.partitions[p])
      continue;

    Model model(parted_msa.model(p));
    assign(model, networkinfo, p);
    parted_msa.model(p, move(model));
  }
}

void assign(Model& model, const NetworkInfo& networkinfo, size_t partition_id)
{
  const pllmod_networkinfo_t& pll_networkinfo = networkinfo.pll_networkinfo();

  if (partition_id >= pll_networkinfo.partition_count)
    throw out_of_range("Partition ID out of range");

  if (!pll_networkinfo.partitions[partition_id])
    return;

  assign(model, pll_networkinfo.partitions[partition_id]);
  model.alpha(pll_networkinfo.alphas[partition_id]);
  if (pll_networkinfo.brlen_scalers)
    model.brlen_scaler(pll_networkinfo.brlen_scalers[partition_id]);
}

void build_clv_network(ProbVector::const_iterator probs, size_t sites, WeightVector::const_iterator weights,
               pll_partition_t* partition, bool normalize, std::vector<double>& clv)
{
  const auto states = partition->states;
  auto clvp = clv.begin();

  for (size_t i = 0; i < sites; ++i)
  {
    if (weights[i] > 0)
    {
      double sum = 0.;
      for (size_t j = 0; j < states; ++j)
        sum += probs[j];

      for (size_t j = 0; j < states; ++j)
      {
        if (sum > 0.)
          clvp[j] =  normalize ? probs[j] / sum : probs[j];
        else
          clvp[j] = 1.0;
      }

      clvp += states;
    }

    /* NB: clv has to be padded, but msa arrays are not! */
    probs += states;
  }

  assert(clvp == clv.end());
}

void set_partition_tips_network(const Options& opts, const MSA& msa, const IDVector& tip_msa_idmap,
                        const PartitionRange& part_region,
                        pll_partition_t* partition, const pll_state_t * charmap)
{
  /* set pattern weights */
  if (!msa.weights().empty())
    pll_set_pattern_weights(partition, msa.weights().data() + part_region.start);

  if (opts.use_prob_msa && msa.probabilistic())
  {
    assert(!(partition->attributes & PLL_ATTRIB_PATTERN_TIP));
    assert(partition->states == msa.states());

    auto normalize = !msa.normalized();
    auto weights_start = msa.weights().cbegin() + part_region.start;

    // we need a libpll function for that!
    auto clv_size = partition->sites * partition->states;
    std::vector<double> tmp_clv(clv_size);
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      auto prob_start = msa.probs(seq_id, part_region.start);
      build_clv_network(prob_start, partition->sites, weights_start, partition, normalize, tmp_clv);
      pll_set_tip_clv(partition, tip_id, tmp_clv.data(), PLL_FALSE);
    }
  }
  else
  {
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      pll_set_tip_states(partition, tip_id, charmap, msa.at(seq_id).c_str() + part_region.start);
    }
  }
}

void set_partition_tips_network(const Options& opts, const MSA& msa, const IDVector& tip_msa_idmap,
                        const PartitionRange& part_region,
                        pll_partition_t* partition, const pll_state_t * charmap,
                        const WeightVector& weights)
{
  assert(!weights.empty());

  const auto pstart = part_region.start;
  const auto pend = part_region.start + part_region.length;

  /* compress weights array by removing all zero entries */
  uintVector comp_weights;
  for (size_t j = pstart; j < pend; ++j)
  {
    if (weights[j] > 0)
      comp_weights.push_back(weights[j]);
  }

  /* now set tip sequences, ignoring all columns with zero weights */
  if (opts.use_prob_msa && msa.probabilistic())
  {
    assert(!(partition->attributes & PLL_ATTRIB_PATTERN_TIP));
    assert(partition->states == msa.states());

    auto normalize = !msa.normalized();
    auto weights_start = msa.weights().cbegin() + part_region.start;

    // we need a libpll function for that!
    auto clv_size = part_region.length * partition->states;
    std::vector<double> tmp_clv(clv_size);
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      auto prob_start = msa.probs(seq_id, part_region.start);
      build_clv_network(prob_start, part_region.length, weights_start, partition, normalize, tmp_clv);
      pll_set_tip_clv(partition, tip_id, tmp_clv.data(), PLL_FALSE);
    }
  }
  else
  {
    std::vector<char> bs_seq(part_region.length);
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      const char * full_seq = msa.at(seq_id).c_str();
      size_t pos = 0;
      for (size_t j = pstart; j < pend; ++j)
      {
        if (weights[j] > 0)
          bs_seq[pos++] = full_seq[j];
      }
      assert(pos == comp_weights.size());

      pll_set_tip_states(partition, tip_id, charmap, bs_seq.data());
    }
  }

  pll_set_pattern_weights(partition, comp_weights.data());
}

pll_partition_t* create_pll_partition_network(const Options& opts, const PartitionInfo& pinfo,
                                      const IDVector& tip_msa_idmap,
                                      const PartitionRange& part_region, const uintVector& weights)
{
  const MSA& msa = pinfo.msa();
  const Model& model = pinfo.model();

  /* part_length doesn't include columns with zero weight */
  const size_t part_length = weights.empty() ? part_region.length :
                             std::count_if(weights.begin() + part_region.start,
                                           weights.begin() + part_region.start + part_region.length,
                                           [](uintVector::value_type w) -> bool
                                             { return w > 0; }
                                           );

  unsigned int attrs = opts.simd_arch;

  if (opts.use_rate_scalers && model.num_ratecats() > 1)
  {
    attrs |= PLL_ATTRIB_RATE_SCALERS;
  }

  if (opts.use_repeats)
  {
    assert(!(opts.use_prob_msa));
    attrs |= PLL_ATTRIB_SITE_REPEATS;
  }
  else if (opts.use_tip_inner)
  {
    assert(!(opts.use_prob_msa));
    // 1) SSE3 tip-inner kernels are not implemented so far, so generic version will be faster
    // 2) same for state-rich models
    if (opts.simd_arch != PLL_ATTRIB_ARCH_SSE && model.num_states() <= 20)
    {
      // TODO: use proper auto-tuning
      const unsigned long min_len_ti = model.num_states() > 4 ? 40 : 100;
      if ((unsigned long) part_length > min_len_ti)
        attrs |= PLL_ATTRIB_PATTERN_TIP;
    }
  }

  // NOTE: if partition is split among multiple threads, asc. bias correction must be applied only once!
  if (model.ascbias_type() == AscBiasCorrection::lewis ||
      (model.ascbias_type() != AscBiasCorrection::none && part_region.master()))
  {
    attrs |=  PLL_ATTRIB_AB_FLAG;
    attrs |= (unsigned int) model.ascbias_type();
  }

  unsigned int max_inner_nodes = msa.size() + MAX_RETICULATIONS;
  unsigned int max_num_branches = msa.size() * 2 + MAX_RETICULATIONS;

  pll_partition_t * partition = pll_partition_create(
	  msa.size(),         /* number of tip sequences */
	  max_inner_nodes,        /* number of CLV buffers */
      model.num_states(),      /* number of states in the data */
      part_length,             /* number of alignment sites/patterns */
      model.num_submodels(),   /* number of different substitution models (LG4 = 4) */
	  max_num_branches,     /* number of probability matrices */
      model.num_ratecats(),    /* number of (GAMMA) rate categories */
	  max_inner_nodes,        /* number of scaling buffers */
      attrs                    /* list of flags (SSE3/AVX, TIP-INNER special cases etc.) */
  );

  libpll_check_error("ERROR creating pll_partition");
  assert(partition);

  if (part_region.master() && !model.ascbias_weights().empty())
    pll_set_asc_state_weights(partition, model.ascbias_weights().data());

  if (part_length == part_region.length)
    set_partition_tips_network(opts, msa, tip_msa_idmap, part_region, partition, model.charmap());
  else
    set_partition_tips_network(opts, msa, tip_msa_idmap, part_region, partition, model.charmap(), weights);

  assign(partition, model);

  return partition;
}
