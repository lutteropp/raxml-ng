#ifndef RAXML_TREEINFO_HPP_
#define RAXML_TREEINFO_HPP_

#include "common.h"
#include "Tree.hpp"
#include "Options.hpp"
#include "loadbalance/PartitionAssignment.hpp"
#include "SPRRoundParams.hpp"

class TreeInfo
{
public:
  TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign);
  TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
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

private:
  pllmod_treeinfo_t * _pll_treeinfo;
  IDSet _parts_master;
  int _brlen_opt_method;
  double _brlen_min;
  double _brlen_max;
  doubleVector _partition_contributions;

  void init(const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
            const std::vector<uintVector>& site_weights);
};

void assign(PartitionedMSA& parted_msa, const TreeInfo& treeinfo);
void assign(Model& model, const TreeInfo& treeinfo, size_t partition_id);


pll_partition_t* create_pll_partition(const Options& opts, const PartitionInfo& pinfo,
                                      const IDVector& tip_msa_idmap,
                                      const PartitionRange& part_region, const uintVector& weights);

#endif /* RAXML_TREEINFO_HPP_ */
