#ifndef RAXML_BOOTSTRAP_BOOTSTRAPNETWORK_HPP_
#define RAXML_BOOTSTRAP_BOOTSTRAPNETWORK_HPP_

#include "SupportNetwork.hpp"

class BootstrapNetwork : public SupportNetwork
{
public:
  BootstrapNetwork (const Network& network);

  virtual
  ~BootstrapNetwork ();

protected:
  virtual void add_network(const pll_unetwork_node_t& root);
};

#endif /* RAXML_BOOTSTRAP_BOOTSTRAPNETWORK_HPP_ */
