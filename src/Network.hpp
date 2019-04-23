#ifndef RAXML_NETWORK_HPP_
#define RAXML_NETWORK_HPP_

#include "common.h"
#include "PartitionedMSA.hpp"

static const int MAX_RETICULATIONS = 64;

// seems to be the only way to have custom deleter for unique_ptr
// without having to specify it every time during object creation
namespace std
{
  template<>
  struct default_delete<pll_rnetwork_t> {
    void operator()(pll_rnetwork_t* ptr) { pll_rnetwork_destroy(ptr, nullptr); }
  };
}

struct NetworkBranch
{
  NetworkBranch() : left_node_id(0), right_node_id(0), length(0.), prob(1.0) {};
  NetworkBranch(size_t left_node_id, size_t right_node_id, double length, double prob) :
    left_node_id(left_node_id), right_node_id(right_node_id), length(length), prob(prob) {};

  size_t left_node_id;
  size_t right_node_id;
  double length;
  double prob;
};

struct NetworkTopology
{
  typedef std::vector<NetworkBranch> edge_container;
  typedef typename edge_container::iterator        iterator;
  typedef typename edge_container::const_iterator  const_iterator;

  //Iterator Compatibility
  iterator begin() { return edges.begin(); }
  iterator end() { return edges.end(); }
  const_iterator begin() const { return edges.cbegin(); }
  const_iterator end() const { return edges.cend(); }
  const_iterator cbegin() { return edges.cbegin(); }
  const_iterator cend() { return edges.cend(); }

  edge_container edges;
  std::vector<doubleVector> brlens;
};

typedef std::unique_ptr<pll_rnetwork_t> PllRNetworkUniquePtr;
typedef std::vector<pll_rnetwork_node_t*> PllNetworkNodeVector;

class BasicNetwork
{
public:
  BasicNetwork(size_t num_tips) : _num_tips(num_tips) {}
  virtual ~BasicNetwork() {}

  bool empty() const { return _num_tips == 0; };
  virtual bool binary() const { return true; };
  virtual size_t num_tips() const { return _num_tips; };
  virtual size_t num_inner_tree() const {return _num_tips - 2;};
  virtual size_t num_reticulations() const {return _num_reticulations;}
  virtual size_t num_inner() const { return num_inner_tree() + num_reticulations(); };
  virtual size_t num_nodes() const { return num_tips() + num_inner(); };
  virtual size_t num_branches() const { return _num_tips + _num_tips - 3 + num_reticulations(); };
  virtual size_t num_tree_branches() const { return _num_tips + _num_tips - 3; };
  virtual size_t num_splits() const { return num_branches() - _num_tips; };

protected:
  size_t _num_tips;
  size_t _num_reticulations;
};

class Network : public BasicNetwork
{
public:
  Network() : BasicNetwork(0), _pll_rnetwork(nullptr) {}
  Network(unsigned int tip_count, const pll_rnetwork_node_t& root) :
    BasicNetwork(tip_count),
    _pll_rnetwork(pll_rnetwork_wrapnetwork(pll_rnetwork_graph_clone(&root))) {}
  Network(const pll_rnetwork_t& pll_rnetwork) :
    BasicNetwork(pll_rnetwork.tip_count), _pll_rnetwork(pll_rnetwork_clone(&pll_rnetwork)) {}
  Network(std::unique_ptr<pll_rnetwork_t>&  pll_rnetwork) :
    BasicNetwork(pll_rnetwork ? pll_rnetwork->tip_count : 0), _pll_rnetwork(pll_rnetwork.release()) {}
  Network(std::unique_ptr<pll_rnetwork_t>&&  pll_rnetwork) :
    BasicNetwork(pll_rnetwork ? pll_rnetwork->tip_count : 0), _pll_rnetwork(pll_rnetwork.release()) {}

  Network (const Network& other);
  Network& operator=(const Network& other);
  Network (Network&& other);
  Network& operator=(Network&& other);

  virtual ~Network();

  static Network buildRandom(size_t num_tips, const char * const* tip_labels, unsigned int random_seed);
  static Network buildRandom(const NameList& taxon_names, unsigned int random_seed);
  static Network buildRandomConstrained(const NameList& taxon_names, unsigned int random_seed,
                                     const Network& constrained_network);
  static Network buildParsimony(const PartitionedMSA& parted_msa, unsigned int random_seed,
                             unsigned int attributes, unsigned int * score = nullptr);
  static Network loadFromFile(const std::string& file_name);

  std::vector<const char*> tip_labels_cstr() const;
  NameList tip_labels_list() const;
  IdNameVector tip_labels() const;
  NameIdMap tip_ids() const;

  NetworkTopology topology() const;
  void topology(const NetworkTopology& topol);

  const std::vector<doubleVector>& partition_brlens() const { return _partition_brlens; }
  const doubleVector& partition_brlens(size_t partition_idx) const;
  void partition_brlens(size_t partition_idx, const doubleVector& brlens);
  void partition_brlens(size_t partition_idx, doubleVector&& brlens);
  void add_partition_brlens(doubleVector&& brlens);

  // TODO: use move semantics to transfer ownership?
  const pll_rnetwork_t& pll_rnetwork() const { return *_pll_rnetwork; }
  pll_rnetwork_t * pll_rnetwork_copy() const;
  void pll_rnetwork(const pll_rnetwork_t&);
  void pll_rnetwork(unsigned int tip_count, const pll_rnetwork_node_t& root);

  const pll_rnetwork_node_t& pll_rnetwork_root() const { return *_pll_rnetwork->root; }
  bool empty() const { return _num_tips == 0; }

  void fix_missing_brlens(double new_brlen = RAXML_BRLEN_DEFAULT);
  void reset_brlens(double new_brlen = RAXML_BRLEN_DEFAULT);
  void apply_partition_brlens(size_t partition_idx);
  void apply_avg_brlens(const doubleVector& partition_contributions);

  void reset_tip_ids(const NameIdMap& label_id_map);
  void reroot(const NameList& outgroup_taxa, bool add_root_node = false);
  void insert_tips_random(const NameList& tip_names, unsigned int random_seed = 0);

public:
  bool binary() const;
  size_t num_inner() const;
  size_t num_inner_tree() const;
  size_t num_reticulations() const;
  size_t num_branches() const;
  size_t num_tree_branches() const;

protected:
  PllRNetworkUniquePtr _pll_rnetwork;
  std::vector<doubleVector> _partition_brlens;

  mutable PllNetworkNodeVector _pll_rnetwork_tips;
  mutable PllNetworkNodeVector _pll_rnetwork_reticulations;
  mutable PllNetworkNodeVector _pll_rnetwork_nodes;

  PllNetworkNodeVector const& tip_nodes() const;
  PllNetworkNodeVector const& reticulation_nodes() const;
  PllNetworkNodeVector const& nodes() const;
};

typedef std::vector<Network> NetworkList;

typedef std::pair<double, NetworkTopology> ScoredNetworkTopology;

class NetworkCollection
{
public:
  typedef std::vector<ScoredNetworkTopology> container_type;
  typedef container_type::const_iterator const_iterator;
  typedef container_type::value_type value_type;

  size_t size() const { return  _networks.size(); }
  bool empty() const { return  _networks.empty(); }
  const_iterator best() const;
  value_type::first_type best_score() const { return best()->first; }
  const value_type::second_type& best_topology() const { return best()->second; }

  const_iterator begin() const { return _networks.cbegin(); }
  const_iterator end() const { return _networks.cend(); }

  void clear() { _networks.clear(); };
  void push_back(double score, const Network& network);
  void push_back(double score, NetworkTopology&& topol);

private:
  container_type _networks;
};

#endif /* RAXML_NETWORK_HPP_ */
