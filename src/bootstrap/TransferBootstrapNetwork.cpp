#include "TransferBootstrapNetwork.hpp"

TransferBootstrapNetwork::TransferBootstrapNetwork(const Network& network, bool naive) :
   SupportNetwork (network), _split_info(nullptr), _naive_method(naive)
{
  assert(num_splits() > 0);
  _node_split_map.resize(num_splits());

  /* extract reference tree splits and add them into hashtable */
  add_network(pll_unetwork_root());

  if (!_naive_method)
  {
    _split_info = pllmod_unetwork_tbe_nature_init((pll_unetwork_t*) &pll_unetwork(), _num_tips,
                                              (const pll_unetwork_node_t**) _node_split_map.data());
  }
}

TransferBootstrapNetwork::~TransferBootstrapNetwork()
{
  if (_split_info)
    free(_split_info);
}

void TransferBootstrapNetwork::add_network(const pll_unetwork_node_t& root)
{
  bool ref_network = (_num_bs_networks == 0);
  doubleVector support(num_splits(), 0.);

  if (ref_network)
  {
    _ref_splits = extract_splits_from_network(root, _node_split_map.data());

    add_splits_to_hashtable(_ref_splits, support, 0);
  }
  else
  {
    assert(_ref_splits);

    auto splits = extract_splits_from_network(root, nullptr);

    // compute TBE
    if (_naive_method)
      pllmod_utree_tbe_naive(_ref_splits.get(), splits.get(), _num_tips, support.data());
    else
    {
      pllmod_unetwork_tbe_nature(_ref_splits.get(), splits.get(), (pll_unetwork_node_t*) &root,
                                               _num_tips, support.data(), _split_info);
    }

    add_splits_to_hashtable(_ref_splits, support, 1);
  }
}
