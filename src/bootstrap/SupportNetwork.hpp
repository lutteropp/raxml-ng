/*
 * SupportNetwork.hpp
 *
 *  Created on: Apr 23, 2019
 *      Author: sarah
 */

#ifndef SRC_BOOTSTRAP_SUPPORTNETWORK_HPP_
#define SRC_BOOTSTRAP_SUPPORTNETWORK_HPP_

#include "../Network.hpp"

typedef std::shared_ptr<pll_split_t> PllSplitSharedPtr;


class SupportNetwork : public Network
{
public:
  SupportNetwork (const Network& network = Network());

  virtual
  ~SupportNetwork ();

  virtual void add_replicate_network(const Network& network);

  void draw_support(bool support_in_pct = true);

protected:
  PllSplitSharedPtr extract_splits_from_network(const pll_unetwork_node_t& root,
                                             pll_unetwork_node_t ** node_split_map);
  void add_splits_to_hashtable(const PllSplitSharedPtr splits,
                               const doubleVector& support, bool update_only);
  void add_network(const Network& network);
  virtual void add_network(const pll_unetwork_node_t& root) = 0;


  void normalize_support_in_hashtable();
  void collect_support();
  virtual bool compute_support();

protected:
  size_t _num_bs_networks;
  bitv_hashtable_t* _pll_splits_hash;
  std::vector<pll_unetwork_node_t*> _node_split_map;
  doubleVector _support;
};

#endif /* SRC_BOOTSTRAP_SUPPORTNETWORK_HPP_ */
