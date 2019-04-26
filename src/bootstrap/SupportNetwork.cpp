/*
 * SupportNetwork.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: sarah
 */

#include "SupportNetwork.hpp"

using namespace std;

char * support_fmt_pct_network(double support)
{
  char *str;
  int size_alloced = asprintf(&str, "%u", (unsigned int) round(support * 100.));

  return size_alloced >= 0 ? str : NULL;
}

char * support_fmt_prop_network(double support)
{
  const unsigned int precision = logger().precision(LogElement::brlen);

  char * str;
  int size_alloced = asprintf(&str, "%.*lf", precision, support);

  return size_alloced >= 0 ? str : NULL;
}


SupportNetwork::SupportNetwork(const Network& network) : Network(network), _num_bs_networks(0)
{
  _pll_splits_hash = nullptr;
}

SupportNetwork::~SupportNetwork ()
{
  if (_pll_splits_hash)
    pllmod_utree_split_hashtable_destroy(_pll_splits_hash);
}

PllSplitSharedPtr SupportNetwork::extract_splits_from_network(const pll_unetwork_node_t& root,
                                                        pll_unetwork_node_t ** node_split_map)
{
  PllSplitSharedPtr splits(pllmod_unetwork_split_create((pll_unetwork_node_t*) &root,
                                                       _num_tips,
													   _num_reticulations,
                                                       node_split_map),
                           pllmod_utree_split_destroy);

  return splits;
}

void SupportNetwork::add_splits_to_hashtable(const PllSplitSharedPtr splits,
                                          const doubleVector& support, bool update_only)
{
  _pll_splits_hash = pllmod_utree_split_hashtable_insert(_pll_splits_hash,
                                                         splits.get(),
                                                         _num_tips,
                                                         num_splits(),
                                                         support.empty() ? nullptr: support.data(),
                                                         update_only);
}


void SupportNetwork::add_network(const Network& network)
{
  add_network(network.pll_unetwork_root());
}

void SupportNetwork::add_replicate_network(const Network& network)
{
  if (network.num_tips() != _num_tips)
    throw runtime_error("Incompatible network!");

  _num_bs_networks++;

  /* extract replicate tree splits and add them into hashtable */
  add_network(network);
//  LOG_DEBUG_TS << "Added replicate trees: " << _num_bs_trees << endl;
}

void SupportNetwork::normalize_support_in_hashtable()
{
  for (unsigned int i = 0; i < _pll_splits_hash->table_size; ++i)
  {
    bitv_hash_entry_t * e =  _pll_splits_hash->table[i];
    while (e != NULL)
    {
      e->support /= _num_bs_networks;
      e = e->next;
    }
  }
}

void SupportNetwork::collect_support()
{
  _support.resize(_pll_splits_hash->entry_count);

  for (unsigned int i = 0; i < _pll_splits_hash->table_size; ++i)
  {
    bitv_hash_entry_t * e =  _pll_splits_hash->table[i];
    while (e != NULL)
    {
      _support[e->bip_number] = e->support;
      e = e->next;
    }
  }
}

bool SupportNetwork::compute_support()
{
  normalize_support_in_hashtable();

  collect_support();

  return true;
}

void SupportNetwork::draw_support(bool support_in_pct)
{
  if (!compute_support())
    return;

  LOG_DEBUG_TS << "Drawing support values on consensus tree..." << endl;

//  printf("\n\n");
//  for (size_t i = 0; i < num_splits(); ++i)
//    printf("node_id %d, split_id %d\n", _node_split_map[i]->node_index, i);
//  printf("\n\n");

  pll_unetwork_node_t ** node_map = _node_split_map.empty() ? nullptr : _node_split_map.data();
  pllmod_unetwork_draw_support(_pll_unetwork.get(), _support.data(), node_map,
                            support_in_pct ? support_fmt_pct_network : support_fmt_prop_network);

  LOG_DEBUG_TS << "Done!" << endl << endl;
}
