#ifndef RAXML_BOOTSTRAP_TRANSFERBOOTSTRAPNETWORK_HPP_
#define RAXML_BOOTSTRAP_TRANSFERBOOTSTRAPNETWORK_HPP_

#include "SupportNetwork.hpp"

class TransferBootstrapNetwork : public SupportNetwork
{
public:
	TransferBootstrapNetwork(const Network& network, bool naive = false);
  virtual ~TransferBootstrapNetwork();

protected:
  virtual void add_network(const pll_unetwork_node_t& root);

protected:
  PllSplitSharedPtr _ref_splits;

private:
  pllmod_tbe_split_info_t * _split_info;
  bool _naive_method;
};

#endif /* RAXML_BOOTSTRAP_TRANSFERBOOTSTRAPNETWORK_HPP_ */
