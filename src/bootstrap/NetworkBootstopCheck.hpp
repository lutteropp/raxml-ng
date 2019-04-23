#ifndef RAXML_BOOTSTRAP_NETWORK_BOOTSTOPCHECK_HPP_
#define RAXML_BOOTSTRAP_NETWORK__BOOTSTOPCHECK_HPP_

#include <bitset>
#include "../Network.hpp"

typedef std::vector<bool> bitVector;
typedef std::vector<bitv_hash_entry_t *> splitEntryVector;

//typedef std::bitset<1000> bitVector;

//class bitVector
//{
//public:
//  bitVector(size_t size): _size(size) { _storage.resize(size / 64, 0); }
//
//  bool operator[](size_t index) const { return _storage[index >> 4] & (1 << (index % 64)); }
//  void set(size_t index) { _storage[index >> 4]  = _storage[index >> 4] | (1 << (index % 64)); }
//
//private:
//  size_t _size;
//  std::vector<size_t> _storage;
//};

class NetworkBootstopCheck
{
protected:
	NetworkBootstopCheck(size_t max_bs_networks);
  virtual ~NetworkBootstopCheck ();

public:
  void add_bootstrap_network(const Network& network);

  bool converged(unsigned long random_seed = 0);

  size_t num_bs_networks() const { return _num_bs_networks; }
  size_t max_bs_networks() const { return _max_bs_networks; }
  void max_bs_networks(size_t val) { if (!_num_bs_networks) _max_bs_networks = val; }

protected:
  size_t _num_bs_networks;
  size_t _max_bs_networks;
  bitv_hashtable_t * _pll_splits_hash;
  std::vector<bitVector> _split_occurence;

  splitEntryVector all_splits();

  virtual bool check_convergence(RandomGenerator& gen) = 0;
};

class NetworkBootstopCheckMRE: public NetworkBootstopCheck
{
public:
	NetworkBootstopCheckMRE(size_t max_bs_networks, double cutoff, size_t num_permutations);
  virtual ~NetworkBootstopCheckMRE ();

  double avg_wrf() const { return _avg_wrf; }
  double avg_pct() const { return _avg_pct; }
  size_t num_better() const { return _num_better; }

protected:
  void mre(splitEntryVector& splits_all, const uintVector& support, splitEntryVector& splits_cons);
  double consensus_wrf_distance(const splitEntryVector& splits1, const splitEntryVector& splits2,
                                const uintVector& support1, const uintVector& support2);

  virtual bool check_convergence(RandomGenerator& gen);

private:
  double _wrf_cutoff;
  size_t _num_permutations;

  double _avg_wrf;
  double _avg_pct;
  size_t _num_better;
};


#endif /* RAXML_BOOTSTRAP_NETWORK_BOOTSTOPCHECK_HPP_ */
