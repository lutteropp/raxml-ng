#include "BootstrapNetwork.hpp"

#include "../common.h"

BootstrapNetwork::BootstrapNetwork (const Network& network) : SupportNetwork(network)
{
  assert(num_splits() > 0);
  _node_split_map.resize(num_splits());

  /* extract reference tree splits and add them into hashtable */
  add_network(pll_unetwork_root());
}

BootstrapNetwork::~BootstrapNetwork ()
{
}

void BootstrapNetwork::add_network(const pll_unetwork_node_t& root)
{
  bool ref_network = (_num_bs_networks == 0);
  pll_unetwork_node_t ** node_split_map = ref_network ? _node_split_map.data() : nullptr;
  int update_only = ref_network ? 0 : 1;
  doubleVector support(num_splits(), ref_network ? 0. : 1.);

  auto splits = extract_splits_from_network(root, node_split_map);

  add_splits_to_hashtable(splits, support, update_only);
}
