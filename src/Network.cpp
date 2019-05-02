#include <algorithm>

#include "Network.hpp"
#include "io/file_io.hpp"

using namespace std;

Network::Network(const Network& other) :
		BasicNetwork(other._num_tips), _pll_unetwork(other.pll_unetwork_copy()), _partition_brlens(other._partition_brlens) {
}

Network::Network(Network&& other) :
		BasicNetwork(other._num_tips), _pll_unetwork(other._pll_unetwork.release()) {
	other._num_tips = 0;
	swap(_pll_unetwork_tips, other._pll_unetwork_tips);
	swap(_pll_unetwork_reticulations, other._pll_unetwork_reticulations);
	swap(_pll_unetwork_nodes, other._pll_unetwork_nodes);
	swap(_partition_brlens, other._partition_brlens);
}

Network& Network::operator=(const Network& other) {
	if (this != &other) {
		_pll_unetwork.reset(other.pll_unetwork_copy());
		_num_tips = other._num_tips;
		_pll_unetwork_tips.clear();
		_pll_unetwork_reticulations.clear();
		_pll_unetwork_nodes.clear();
		_partition_brlens = other._partition_brlens;
	}

	return *this;
}

Network& Network::operator=(Network&& other) {
	if (this != &other) {
		_num_tips = 0;
		_pll_unetwork_tips.clear();
		_pll_unetwork_reticulations.clear();
		_pll_unetwork_nodes.clear();
		_pll_unetwork.release();
		_partition_brlens.clear();

		swap(_num_tips, other._num_tips);
		swap(_pll_unetwork, other._pll_unetwork);
		swap(_pll_unetwork_tips, other._pll_unetwork_tips);
		swap(_pll_unetwork_reticulations, other._pll_unetwork_reticulations);
		swap(_pll_unetwork_nodes, other._pll_unetwork_nodes);
		swap(_partition_brlens, other._partition_brlens);
	}

	return *this;
}

Network::~Network() {
}

size_t Network::num_inner_induced_tree() const {
	return _pll_unetwork ? _pll_unetwork->inner_tree_count : BasicNetwork::num_inner_induced_tree();
}

size_t Network::num_inner() const {
	return _pll_unetwork ? _pll_unetwork->inner_tree_count + _pll_unetwork->reticulation_count : BasicNetwork::num_inner();
}

size_t Network::num_reticulations() const {
	return _pll_unetwork ? _pll_unetwork->reticulation_count : BasicNetwork::num_reticulations();
}

size_t Network::num_branches() const {
	return _pll_unetwork ? _pll_unetwork->edge_count : BasicNetwork::num_branches();
}

size_t Network::num_branches_induced_tree() const {
	return _pll_unetwork ? _pll_unetwork->tree_edge_count : BasicNetwork::num_branches_induced_tree();
}

bool Network::binary() const {
	return _pll_unetwork ? _pll_unetwork->binary : BasicNetwork::binary();
}

pll_unetwork_t * Network::pll_unetwork_copy() const {
	return _pll_unetwork ? pll_unetwork_clone(_pll_unetwork.get()) : nullptr;
}

void Network::pll_unetwork(const pll_unetwork_t& network) {
	_num_tips = network.tip_count;
	_pll_unetwork.reset(pll_unetwork_clone(&network));
	_pll_unetwork_tips.clear();
	_pll_unetwork_reticulations.clear();
	_pll_unetwork_nodes.clear();
}

Network Network::buildRandom(size_t num_tips, const char * const * tip_labels, unsigned int random_seed) {
	PllUNetworkUniquePtr pll_unetwork(pllmod_unetwork_create_random(num_tips, tip_labels, random_seed));

	libpll_check_error("ERROR building random network");
	assert(pll_unetwork);

	return Network(pll_unetwork);
}

Network Network::buildRandom(const NameList& taxon_names, unsigned int random_seed) {
	std::vector<const char*> tip_labels(taxon_names.size(), nullptr);
	for (size_t i = 0; i < taxon_names.size(); ++i)
		tip_labels[i] = taxon_names[i].data();

	return Network::buildRandom(taxon_names.size(), (const char * const *) tip_labels.data(), random_seed);
}

Network Network::buildRandomConstrained(const NameList& taxon_names, unsigned int random_seed, const Network& constrained_network) {
	PllUNetworkUniquePtr pll_unetwork(pllmod_unetwork_resolve_multi(&constrained_network.pll_unetwork(), random_seed, nullptr));

	if (!pll_unetwork) {
		assert (pll_errno);
		libpll_check_error("ERROR in building a randomized constrained network");
	}

	Network network(pll_unetwork);

	if (taxon_names.size() > network.num_tips()) {
		// constraint network is not comprehensive -> add free taxa
		auto free_tip_count = taxon_names.size() - network.num_tips();
		auto cons_tips = network.tip_ids();
		NameList free_tips;
		free_tips.reserve(free_tip_count);

		for (const auto& t : taxon_names) {
			if (!cons_tips.count(t))
				free_tips.push_back(t);
		}

		assert(free_tips.size() == free_tip_count);

		network.insert_tips_random(free_tips, random_seed);
	}

//  pll_rnetwork_check_integrity(&network.pll_rnetwork());

	return network;
}

Network Network::buildParsimony(const PartitionedMSA& parted_msa, unsigned int random_seed, unsigned int attributes, unsigned int * score) {
	unsigned int lscore;
	unsigned int *pscore = score ? score : &lscore;

	LOG_DEBUG << "Parsimony seed: " << random_seed << endl;

	auto taxon_names = parted_msa.taxon_names();
	auto taxon_count = taxon_names.size();
	std::vector<const char*> tip_labels(taxon_names.size(), nullptr);
	for (size_t i = 0; i < taxon_names.size(); ++i)
		tip_labels[i] = taxon_names[i].data();

	// create pll_partitions
	std::vector<pll_partition*> pars_partitions(parted_msa.part_count(), nullptr);
	size_t i = 0;
	for (const auto& pinfo : parted_msa.part_list()) {
		const auto& model = pinfo.model();
		const auto& msa = pinfo.msa();

		auto partition = pll_partition_create(taxon_count, 0, /* number of CLVs */
		model.num_states(), msa.length(), 1, 1, /* pmatrix count */
		1, /* rate_cats */
		0, /* scale buffers */
		attributes);

		/* set pattern weights */
		if (!msa.weights().empty())
			pll_set_pattern_weights(partition, msa.weights().data());

		/* set tip states */
		for (size_t j = 0; j < msa.size(); ++j) {
			pll_set_tip_states(partition, j, model.charmap(), msa.at(j).c_str());
		}

		pars_partitions[i++] = partition;
	}
	assert(i == pars_partitions.size());

	PllUNetworkUniquePtr pll_unetwork(
			pllmod_unetwork_create_parsimony_multipart(taxon_count, (char* const *) tip_labels.data(), pars_partitions.size(),
					pars_partitions.data(), random_seed, pscore));

	for (auto p : pars_partitions)
		pll_partition_destroy(p);

	libpll_check_error("ERROR building parsimony network");
	assert(pll_unetwork);

	return Network(pll_unetwork);
}

Network Network::loadFromFile(const std::string& file_name) {
	Network network;
	NetworkNewickStream ns(file_name, std::ios::in);

	ns >> network;

	return network;
}

PllNetworkNodeVector const& Network::tip_nodes() const {
	if (_pll_unetwork_tips.empty() && _num_tips > 0) {
		assert(_num_tips == _pll_unetwork->tip_count);

		_pll_unetwork_tips.assign(_pll_unetwork->nodes, _pll_unetwork->nodes + _pll_unetwork->tip_count);
	}

	return _pll_unetwork_tips;
}

PllNetworkNodeVector const& Network::nodes() const {
	if (_pll_unetwork_nodes.empty() && _num_tips > 0) {
		assert(_num_tips == _pll_unetwork->tip_count);

		_pll_unetwork_nodes.assign(_pll_unetwork->nodes,
				_pll_unetwork->nodes + _pll_unetwork->tip_count + _pll_unetwork->inner_tree_count + _pll_unetwork->reticulation_count);
	}

	return _pll_unetwork_tips;
}

PllNetworkNodeVector const& Network::reticulation_nodes() const {
	if (_pll_unetwork_reticulations.empty() && _num_reticulations > 0) {
		assert(_num_reticulations == _pll_unetwork->reticulation_count);

		// TODO: _pll_rnetwork_reticulations.assign(...);
	}

	return _pll_unetwork_reticulations;
}

std::vector<const char*> Network::tip_labels_cstr() const {
	std::vector<const char*> result;

	if (!empty()) {
		result.resize(_num_tips, nullptr);
		for (auto const& node : tip_nodes())
			result[node->clv_index] = node->label;
	}

	return result;
}

NameList Network::tip_labels_list() const {
	NameList result;

	if (!empty()) {
		result.resize(_num_tips);
		for (auto const& node : tip_nodes())
			result[node->clv_index] = string(node->label);
	}

	return result;
}

IdNameVector Network::tip_labels() const {
	IdNameVector result;
	for (auto const& node : tip_nodes())
		result.emplace_back(node->clv_index, string(node->label));

	assert(!result.empty());

	return result;
}

NameIdMap Network::tip_ids() const {
	NameIdMap result;
	for (auto const& node : tip_nodes())
		result.emplace(string(node->label), node->clv_index);

	assert(!result.empty());

	return result;
}

void Network::insert_tips_random(const NameList& tip_names, unsigned int random_seed) {
	_pll_unetwork_tips.clear();

	std::vector<const char*> tip_labels(tip_names.size(), nullptr);
	for (size_t i = 0; i < tip_names.size(); ++i)
		tip_labels[i] = tip_names[i].data();

	int retval = pllmod_unetwork_extend_random(_pll_unetwork.get(), tip_labels.size(), (const char * const *) tip_labels.data(),
			random_seed);

	if (retval)
		_num_tips = _pll_unetwork->tip_count;
	else {
		assert (pll_errno);
		libpll_check_error("ERROR in randomized network extension");
	}
}

void Network::reset_tip_ids(const NameIdMap& label_id_map) {
	if (label_id_map.size() < _num_tips)
		throw invalid_argument("Invalid map size");

	for (auto& node : tip_nodes()) {
		const unsigned int tip_id = label_id_map.at(node->label);
		node->clv_index = node->node_index = tip_id;
	}
}

void Network::fix_missing_brlens(double new_brlen) {
	pllmod_unetwork_set_length_recursive(_pll_unetwork.get(), new_brlen, 1);
}

void Network::reset_brlens(double new_brlen) {
	pllmod_unetwork_set_length_recursive(_pll_unetwork.get(), new_brlen, 0);
}

NetworkTopology Network::topology() const {
	NetworkTopology topol;

	  topol.edges.resize(num_branches());

	  size_t branches = 0;
	  for (auto n: subnodes())
	  {
	    if (n->node_index < n->back->node_index /*&& n->active*/)
	    {
	      assert(n->pmatrix_index < topol.edges.size());
	      topol.edges.at(n->pmatrix_index) = NetworkBranch(n->node_index, n->back->node_index, n->length, n->prob);
	      branches++;
	    }
	  }
	  topol.brlens = _partition_brlens;

	//  for (auto& branch: topol.edges)
	//    printf("%u %u %lf\n", branch.left_node_id, branch.right_node_id, branch.length);

	  assert(branches == num_branches());

	  return topol;
}

void Network::topology(const NetworkTopology& topol) {
	if (topol.edges.size() != num_branches())
		throw runtime_error("Incompatible topology!");

	auto allnodes = nodes();
	unsigned int pmatrix_index = 0;
	for (const auto& branch : topol) {
		pll_unetwork_node_t * left_node = allnodes.at(branch.left_node_id);
		pll_unetwork_node_t * right_node = allnodes.at(branch.right_node_id);
		pllmod_unetwork_connect_nodes(left_node, right_node, branch.length, branch.prob);

		// important: make sure all branches have distinct pmatrix indices!
		left_node->pmatrix_index = right_node->pmatrix_index = pmatrix_index;

		pmatrix_index++;
//    printf("%u %u %lf %d  (%u - %u) \n", branch.left_node_id, branch.right_node_id,
//           branch.length, left_node->pmatrix_index, left_node->clv_index, right_node->clv_index);
	}

	_partition_brlens = topol.brlens;

	assert(pmatrix_index == num_branches());
}

const doubleVector& Network::partition_brlens(size_t partition_idx) const {
	return _partition_brlens.at(partition_idx);
}

void Network::partition_brlens(size_t partition_idx, const doubleVector& brlens) {
	_partition_brlens.at(partition_idx) = brlens;
}

void Network::partition_brlens(size_t partition_idx, doubleVector&& brlens) {
	_partition_brlens.at(partition_idx) = brlens;
}

void Network::add_partition_brlens(doubleVector&& brlens) {
	_partition_brlens.push_back(brlens);
}

void Network::apply_partition_brlens(size_t partition_idx) {
	if (partition_idx >= _partition_brlens.size())
		throw out_of_range("Partition ID out of range");

	const auto brlens = _partition_brlens.at(partition_idx);
	for (auto n : nodes()) {
		n->length = brlens[n->pmatrix_index];
	}
}

void Network::apply_avg_brlens(const doubleVector& partition_contributions) {
	assert(!_partition_brlens.empty() && partition_contributions.size() == _partition_brlens.size());

	const auto allnodes = nodes();

	for (auto n : allnodes)
		n->length = 0;

	for (size_t p = 0; p < _partition_brlens.size(); ++p) {
		const auto brlens = _partition_brlens[p];
		const auto w = partition_contributions[p];
		for (auto n : allnodes)
			n->length += brlens[n->pmatrix_index] * w;
	}
}

PllNetworkNodeVector Network::subnodes() const
{
	PllNetworkNodeVector subnodes;

  if (_num_tips > 0)
  {
    subnodes.resize(num_subnodes());

    for (size_t i = 0; i < _pll_unetwork->tip_count + _pll_unetwork->inner_tree_count + _pll_unetwork->reticulation_count; ++i)
    {
      auto start = _pll_unetwork->nodes[i];
      auto node = start;
      do
      {
        subnodes[node->node_index] = node;
        node = node->next;
      }
      while (node && node != start);
    }
  }

  return subnodes;
}

void Network::reroot(const NameList& outgroup_taxa, bool add_root_node) {
	// collect tip node indices
	NameIdMap name_id_map;
	for (auto const& node : tip_nodes())
		name_id_map.emplace(string(node->label), node->node_index);

	// find tip ids for outgroup taxa
	uintVector tip_ids;
	for (const auto& label : outgroup_taxa) {
		const auto tip_id = name_id_map.at(label);
		tip_ids.push_back(tip_id);
	}

	// re-root network with the outgroup
	int res = pllmod_unetwork_outgroup_root(_pll_unetwork.get(), tip_ids.data(), tip_ids.size(), add_root_node);

	if (!res)
		libpll_check_error("Unable to reroot network");
}

NetworkCollection::const_iterator NetworkCollection::best() const {
	return std::max_element(_networks.cbegin(), _networks.cend(), [](const value_type& a, const value_type& b) -> bool
	{	return a.first < b.first;});
}

void NetworkCollection::push_back(double score, const Network& network) {
	_networks.emplace_back(score, network.topology());
}

void NetworkCollection::push_back(double score, NetworkTopology&& topol) {
	_networks.emplace_back(score, topol);
}
