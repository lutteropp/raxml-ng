#ifndef RAXML_TREEINFO_HPP_
#define RAXML_TREEINFO_HPP_

#include "common.h"
#include "Tree.hpp"
#include "Options.hpp"
#include "AncestralStates.hpp"
#include "loadbalance/PartitionAssignment.hpp"

struct spr_round_params
{
  bool thorough;
  int radius_min;
  int radius_max;
  int ntopol_keep;
  double subtree_cutoff;
  cutoff_info_t cutoff_info;

  void reset_cutoff_info(double loglh)
  {
    cutoff_info.lh_dec_count = 0;
    cutoff_info.lh_dec_sum = 0.;
    cutoff_info.lh_cutoff = loglh / -1000.0;
  }
};

class TreeInfo
{
public:
  TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign);
  TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
            const std::vector<uintVector>& site_weights);

  TreeInfo (const Options &opts, const std::vector<doubleVector>& partition_brlens, pllmod_treeinfo_t* base_treeinfo,
		      const PartitionedMSA& parted_msa,
              const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign);
  TreeInfo (const Options &opts, const std::vector<doubleVector>& partition_brlens, pllmod_treeinfo_t* base_treeinfo,
		      const PartitionedMSA& parted_msa,
              const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
              const std::vector<uintVector>& site_weights);
  TreeInfo (const Options &opts, pllmod_treeinfo_t* base_treeinfo,
		      const PartitionedMSA& parted_msa,
              const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign);
  TreeInfo (const Options &opts, pllmod_treeinfo_t* base_treeinfo,
		      const PartitionedMSA& parted_msa,
              const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
              const std::vector<uintVector>& site_weights);

  virtual
  ~TreeInfo ();

  const pllmod_treeinfo_t& pll_treeinfo() const { return *_pll_treeinfo; }
  const pll_unode_t& pll_utree_root() const { assert(_pll_treeinfo); return *_pll_treeinfo->root; }

  Tree tree() const;
  Tree tree(size_t partition_id) const;
  void tree(const Tree& tree);

  /* in parallel mode, partition can be share among multiple threads and TreeInfo objects;
   * this method returns list of partition IDs for which this thread is designated as "master"
   * and thus responsible for e.g. sending model parameters to the main thread. */
  const IDSet& parts_master() const { return _parts_master; }

  void model(size_t partition_id, const Model& model);

  void set_topology_constraint(const Tree& cons_tree);

  double loglh(bool incremental = false);
  double optimize_params(int params_to_optimize, double lh_epsilon);
  double optimize_params_all(double lh_epsilon)
  { return optimize_params(PLLMOD_OPT_PARAM_ALL, lh_epsilon); } ;
  double optimize_model(double lh_epsilon)
  { return optimize_params(PLLMOD_OPT_PARAM_ALL & ~PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE, lh_epsilon); } ;
  double optimize_branches(double lh_epsilon, double brlen_smooth_factor);
  double spr_round(spr_round_params& params);
  void compute_ancestral(const AncestralStatesSharedPtr& ancestral,
                         const PartitionAssignment& part_assign);

  double (*opt_brlen_function)(pllmod_treeinfo_t *, // treeinfo
                                          double, // min_brlen
                                          double, // max_brlen
                                          double, // lh_epsilon
                                          int, // max_iters
                                          int, // opt_method
                                          int); // radius
  double (*spr_round_function)(pllmod_treeinfo_t *, // treeinfo
                                            unsigned int, // radius_min
                                            unsigned int, // radius_max
                                            unsigned int, // ntopol_keep
                                            pll_bool_t, // thorough
                                            int, // brlen_opt_method
                                            double, // bl_min
                                            double, // bl_max
                                            int, // smoothings
                                            double, // epsilon
                                            cutoff_info_t *, // cutoff_info
                                            double); // subtree_cutoff
  pllmod_ancestral_t * (*compute_ancestral_function)(pllmod_treeinfo_t *); // treeinfo
private:
  pllmod_treeinfo_t * _pll_treeinfo;
  IDSet _parts_master;
  int _brlen_opt_method;
  double _brlen_min;
  double _brlen_max;
  bool _check_lh_impr;
  doubleVector _partition_contributions;

  pllmod_treeinfo_t* create_base_treeinfo(const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa);

  void init(const Options &opts, const std::vector<doubleVector>& partition_brlens, pllmod_treeinfo_t* base_treeinfo, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
            const std::vector<uintVector>& site_weights);

  void assert_lh_improvement(double old_lh, double new_lh, const std::string& where = "");
};

void assign(PartitionedMSA& parted_msa, const TreeInfo& treeinfo);
void assign(Model& model, const TreeInfo& treeinfo, size_t partition_id);


pll_partition_t* create_pll_partition(const Options& opts, const PartitionInfo& pinfo,
                                      const IDVector& tip_msa_idmap,
                                      const PartitionRange& part_region, const uintVector& weights);

#endif /* RAXML_TREEINFO_HPP_ */
