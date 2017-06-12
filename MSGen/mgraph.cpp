#include "mgraph.h"
#include <algorithm>

MuCluster::MuCluster(MSGraph & g, ID id, const BitSeq & vec, MuClass & _cls)
	: graph(g), cluster_id(id), _class(_cls), score_vector(vec) {
	BitSeq::size_t k = 0, n = score_vector.bit_number();

	score_degree = 0;
	while (k < n) {
		if (score_vector.get_bit(k++) == BIT_1)
			score_degree++;
	}

	in_port = new MuSubsumePort();
	ou_port = new MuSubsumePort();
}
MuCluster::~MuCluster() {
	delete in_port;
	delete ou_port;
}
void MuCluster::link_to(MuCluster & trg) {
	if (&trg == this) {
		CError error(CErrorType::InvalidArguments, "MuCluster::link_to", "Invalid link: self-connect");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		ou_port->link(*this, trg);
		trg.in_port->link(*this, trg);
	}
}
const std::set<MuCluster *> & MuHierarchy::get_clusters_of(size_t deg) const {
	if (degree_map.count(deg) == 0) {
		CError error(CErrorType::InvalidArguments, "MuHierarchy::get_cluster_of", 
			"Undefined degree (" + std::to_string(deg) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = degree_map.find(deg);
		return *(iter->second);
	}
}
const std::set<MuCluster *> & MuHierarchy::get_clusters_at(size_t index) const {
	if (index >= degree_list.size()) {
		CError error(CErrorType::OutOfIndex, 
			"MuHierarchy::get_clusters_at", 
			"Invalid index (" + std::to_string(index) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::OutOfIndex);
	}
	else {
		size_t deg = degree_list[index];
		auto iter = degree_map.find(deg);
		return *(iter->second);
	}
}
void MuHierarchy::add(MuCluster & cluster) {
	/* get the degree of cluster */
	size_t deg = cluster.get_score_degree();
	
	/* get the target set of clusters in hierarchy */
	std::set<MuCluster *> * clusters;
	if (degree_map.count(deg) == 0) {
		clusters = new std::set<MuCluster *>();
		degree_map[deg] = clusters; 
		degree_list.push_back(deg);
	}
	else {
		auto iter = degree_map.find(deg);
		clusters = iter->second;
	}

	/* insert cluster to the level */ 
	clusters->insert(&cluster);
}
void MuHierarchy::clear() {
	auto beg = degree_map.begin();
	auto end = degree_map.end();
	while (beg != end)
		delete (beg++)->second;

	degree_list.clear(); 
	degree_map.clear();
}
void MuHierarchy::sort() {
	std::sort(degree_list.begin(), degree_list.end());
}

MuCluster & MSGraph::get_cluster(MuCluster::ID cid) const {
	if (cid >= clusters.size()) {
		CError error(CErrorType::OutOfIndex, "MSGraph::get_cluster", 
			"Invalid cluster-id (" + std::to_string(cid) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::OutOfIndex);
	}
	else return *(clusters[cid]);
}
void MSGraph::clear() {
	roots.clear();
	leafs.clear();
	hierarchy.clear();
	index.clear();

	for (size_t i = 0; i < clusters.size(); i++)
		delete clusters[i];
	clusters.clear();
}
void MSGraph::add(MuClassSet & class_set) {
	/* initialization */ clear();

	/* create clusters for each class in the set */
	const std::map<MuFeature, MuClass *> &
		classes = class_set.get_classes();
	auto cbeg = classes.begin(), cend = classes.end();
	while (cbeg != cend) {
		/* get next class and create its cluster in the list */
		MuClass & _class = *((cbeg++)->second);
		BitSeq * score_vec = (BitSeq *)(_class.get_feature());
		MuCluster * cluster = new MuCluster(*this, clusters.size(), *score_vec, _class);
		clusters.push_back(cluster); hierarchy.add(*cluster);
	}

	/* update index for graph */
	Mutant::ID mid, mnum = class_set.get_mutants().get_space().number_of_mutants();
	for (mid = 0; mid < mnum; mid++) {
		size_t k = 0, n = clusters.size();
		while (k < n) {
			if (clusters[k]->has_mutant(mid)) {
				if (index.count(mid) == 0) {
					index[mid] = clusters[k];
				}
				else {
					CError error(CErrorType::Runtime, "MSGraph::add", 
						"Duplicated mutants in two classes (" + std::to_string(mid) + ")");
					CErrorConsumer::consume(error); exit(CErrorType::Runtime);
				}
			}
			k++;
		}
	} /* end for */

	/* sort hierarchy */ hierarchy.sort();
}
void MSGraph::update_roots_and_leafs() {
	/* initialization */
	roots.clear(); leafs.clear();

	/* iterate each cluster */
	size_t i, n = clusters.size();
	for (i = 0; i < n; i++) {
		MuCluster & cluster = *(clusters[i]);

		if (cluster.get_in_port().get_degree() == 0)
			roots.insert(&cluster);
		if (cluster.get_ou_port().get_degree() == 0)
			leafs.insert(&cluster);
	}
}



