#include "ConsensusNetwork.hpp"

using namespace std;

ConsensusNetwork::ConsensusNetwork (const NetworkList& networks, unsigned int consense_cutoff) : SupportNetwork()
{
  assert(consense_cutoff <= 100);
  assert(!networks.empty());

  _cutoff = ((double) consense_cutoff) / 100.;
  pll_unetwork(networks[0].pll_unetwork());

  LOG_DEBUG_TS << "Extracting splits from replicate trees..." << endl;

  for (const auto& t: networks)
    add_replicate_network(t);

  assert(_num_bs_networks == networks.size());

  LOG_DEBUG_TS << "Extracted " << _pll_splits_hash->entry_count << " splits from "
               << _num_bs_networks << " networks." <<  endl;
}

ConsensusNetwork::~ConsensusNetwork ()
{
}

void ConsensusNetwork::add_network(const pll_unetwork_node_t& root)
{
  pll_unetwork_node_t ** node_split_map = nullptr;
  doubleVector support;
  int update_only = 0;

  auto splits = extract_splits_from_network(root, node_split_map);

  add_splits_to_hashtable(splits, support, update_only);
}

bool ConsensusNetwork::compute_support()
{
  LOG_DEBUG_TS << "Building consensus split system..." << endl;

  /* normalize support values */
  normalize_support_in_hashtable();

  /* build final split system */
  pll_split_system_t * split_system = pllmod_utree_split_consensus(_pll_splits_hash,
                                                                   _num_tips,
                                                                   _cutoff);

  if (!split_system)
    libpll_check_error("Failed to create consensus tree (pllmod_unetwork_split_consensus)");

  LOG_DEBUG_TS << "Split system size: " << split_system->split_count << endl;

  LOG_DEBUG_TS << "Building consensus tree topology..." << endl;

  /* build network from splits */
  pll_consensus_unetwork_t * cons_network = pllmod_unetwork_from_splits(split_system,
                                                               _num_tips,
                                                               (char * const *) tip_labels_cstr().data());

  if (!cons_network)
    libpll_check_error("Failed to create consensus network (pllmod_unetwork_from_splits)");

  LOG_DEBUG_TS << "Consensus network has " << to_string(cons_network->branch_count)
               << " internal branches." << endl;

  /* set consensus network topology */
  pll_unetwork(_num_tips, *cons_network->network);

  /* map pll_unodes to splits */
  _node_split_map.resize(_pll_unetwork->inner_tree_count + _pll_unetwork->reticulation_count);
  _support.resize(_pll_unetwork->inner_tree_count + _pll_unetwork->reticulation_count);
  for (unsigned int i = 0; i < _pll_unetwork->inner_tree_count + _pll_unetwork->reticulation_count; ++i)
  {
    auto node = _pll_unetwork->nodes[_pll_unetwork->tip_count + i];
    assert(node->data);
    _node_split_map[i] = node;
    _support[i] = ((pll_consensus_data_t *) node->data)->support;
  }

  pllmod_utree_split_system_destroy(split_system);
  pllmod_unetwork_consensus_destroy(cons_network);

  return true;
}


