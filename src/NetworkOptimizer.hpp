#ifndef RAXML_NETWORK_OPTIMIZER_H_
#define RAXML_NETWORK_OPTIMIZER_H_

#include "NetworkInfo.hpp"
#include "NetworkCheckpoint.hpp"

class NetworkOptimizer
{
public:
	NetworkOptimizer (const Options& opts);
  virtual
  ~NetworkOptimizer ();

  double optimize_model(NetworkInfo& networkinfo, double lh_epsilon);
  double optimize_model(NetworkInfo& networkinfo) { return optimize_model(networkinfo, _lh_epsilon); };
  double optimize_topology(NetworkInfo& networkinfo, NetworkCheckpointManager& cm);
  double evaluate(NetworkInfo& treeinfo, NetworkCheckpointManager& cm);
private:
  double _lh_epsilon;
  int _spr_radius;
  double _spr_cutoff;
};

#endif /* RAXML_NETWORK_OPTIMIZER_H_ */
