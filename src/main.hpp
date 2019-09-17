/*
 * main.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: Sarah Lutteropp
 */

#ifndef SRC_MAIN_HPP_
#define SRC_MAIN_HPP_

#include <stddef.h>
#include <map>
#include <memory>
#include <string>

#include "AncestralStates.hpp"
#include "bootstrap/BootstrapGenerator.hpp"
#include "bootstrap/SupportTree.hpp"
#include "Checkpoint.hpp"
#include "loadbalance/PartitionAssignment.hpp"
#include "MSA.hpp"
#include "Options.hpp"
#include "PartitionedMSA.hpp"
#include "PartitionedMSAView.hpp"
#include "types.hpp"
#include "Tree.hpp"

class BootstopCheckMRE;
class ConsensusTree;
class LoadBalancer;
class NewickStream;
class RFDistCalculator;

using namespace std;

struct RaxmlInstance
{
  Options opts;
  shared_ptr<PartitionedMSA> parted_msa;
  unique_ptr<PartitionedMSA> parted_msa_parsimony;
  TreeList start_trees;
  BootstrapReplicateList bs_reps;
  TreeList bs_start_trees;
  PartitionAssignmentList proc_part_assign;
  unique_ptr<LoadBalancer> load_balancer;
  map<BranchSupportMetric, shared_ptr<SupportTree> > support_trees;
  shared_ptr<ConsensusTree> consens_tree;

  // bootstopping convergence test, only autoMRE is supported for now
  unique_ptr<BootstopCheckMRE> bootstop_checker;

  // mapping taxon name -> tip_id/clv_id in the tree
  NameIdMap tip_id_map;

  // mapping tip_id in the tree (array index) -> sequence index in MSA
  IDVector tip_msa_idmap;

 // unique_ptr<TerraceWrapper> terrace_wrapper;

//  unique_ptr<RandomGenerator> starttree_seed_gen;
//  unique_ptr<RandomGenerator> bootstrap_seed_gen;

  unique_ptr<NewickStream> start_tree_stream;

  /* this is just a dummy random tree used for convenience, e,g, if we need tip labels or
   * just 'any' valid tree for the alignment at hand */
  Tree random_tree;

  /* topological constraint */
  Tree constraint_tree;

  unique_ptr<RFDistCalculator> dist_calculator;
  AncestralStatesSharedPtr ancestral_states;
};

void print_banner();
void init_part_info(RaxmlInstance& instance);
void print_reduced_msa(const RaxmlInstance& instance, const PartitionedMSAView& reduced_msa_view);
bool check_msa_global(const MSA& msa);
bool check_msa(RaxmlInstance& instance);
size_t total_free_params(const RaxmlInstance& instance);
void check_models(const RaxmlInstance& instance);
void check_tree(const PartitionedMSA& msa, const Tree& tree);
void check_options(RaxmlInstance& instance);
void load_msa(RaxmlInstance& instance);
void load_parted_msa(RaxmlInstance& instance);
void prepare_tree(const RaxmlInstance& instance, Tree& tree);
Tree generate_tree(const RaxmlInstance& instance, StartingTree type);
void load_start_trees(RaxmlInstance& instance, CheckpointManager& cm);
void load_checkpoint(RaxmlInstance& instance, CheckpointManager& cm);
void load_constraint(RaxmlInstance& instance);
void build_parsimony_msa(RaxmlInstance& instance);
void build_start_trees(RaxmlInstance& instance, size_t skip_trees);
void balance_load(RaxmlInstance& instance);
void balance_load(RaxmlInstance& instance, WeightVectorList part_site_weights);
void generate_bootstraps(RaxmlInstance& instance, const Checkpoint& checkp);
void init_ancestral(RaxmlInstance& instance);
void reroot_tree_with_outgroup(const Options& opts, Tree& tree, bool add_root_node);
void postprocess_tree(const Options& opts, Tree& tree);
void draw_bootstrap_support(RaxmlInstance& instance, Tree& ref_tree, const TreeCollection& bs_trees);
bool check_bootstop(const RaxmlInstance& instance, const TreeCollection& bs_trees,
                    bool print = false);
TreeCollection read_newick_trees(Tree& ref_tree,
                                 const std::string& fname, const std::string& tree_kind);
TreeCollection read_bootstrap_trees(const RaxmlInstance& instance, Tree& ref_tree);
void read_multiple_tree_files(RaxmlInstance& instance);
void command_bootstop(RaxmlInstance& instance);
void command_support(RaxmlInstance& instance);
void command_rfdist(RaxmlInstance& instance);
void command_consense(RaxmlInstance& instance);
void command_bsmsa(RaxmlInstance& instance, const Checkpoint& checkp);
void check_terrace(const RaxmlInstance& instance, const Tree& tree);
void save_ml_trees(const Options& opts, const Checkpoint& checkp);
void print_ic_scores(const RaxmlInstance& instance, double loglh);
void print_final_output(const RaxmlInstance& instance, const Checkpoint& checkp);
void print_resources(const RaxmlInstance& instance);
void thread_main(RaxmlInstance& instance, CheckpointManager& cm);
void master_main(RaxmlInstance& instance, CheckpointManager& cm);
int clean_exit(int retval);
int internal_main(int argc, char** argv, void* comm);



#endif /* SRC_MAIN_HPP_ */
