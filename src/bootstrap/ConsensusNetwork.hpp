/*
 * ConsensusNetwork.hpp
 *
 *  Created on: Apr 23, 2019
 *      Author: sarah
 */

#ifndef SRC_BOOTSTRAP_CONSENSUSNETWORK_HPP_
#define SRC_BOOTSTRAP_CONSENSUSNETWORK_HPP_

#include "SupportNetwork.hpp"

class ConsensusNetwork : public SupportNetwork
{
public:
	ConsensusNetwork (const NetworkList& networks, unsigned int consense_cutoff);
  virtual
  ~ConsensusNetwork ();

protected:
  virtual void add_network(const pll_rnetwork_node_t& root);
  virtual bool compute_support();

private:
   double _cutoff;
};


#endif /* SRC_BOOTSTRAP_CONSENSUSNETWORK_HPP_ */
