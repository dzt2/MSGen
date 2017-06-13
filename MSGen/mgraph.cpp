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
MuCluster & MSGraph::get_cluster_of(Mutant::ID mid) const {
	if (index.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"MSGraph::get_cluster_of", 
			"Invalid mutant (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error); 
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = index.find(mid);
		return *(iter->second);
	}
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
void MSGraph::clear_edges() {
	size_t i, n = clusters.size();
	for (i = 0; i < n; i++) {
		MuCluster & cluster = *(clusters[i]);
		cluster.get_in_port().clear();
		cluster.get_ou_port().clear();
	}
}
void MSGraph::build(MuClassSet & class_set) {
	/* initialization */ 
	clear(); _class_set = &class_set;

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

_MSG_VSpace_down_top::~_MSG_VSpace_down_top() {
	while (!vqueue.empty()) vqueue.pop();
	while (!vcache.empty()) vcache.pop();
	visitset.clear();
}
void _MSG_VSpace_down_top::initial() {
	while (!vqueue.empty()) vqueue.pop();
	while (!vcache.empty()) vcache.pop();
	visitset.clear();

	auto beg = leafs.begin();
	auto end = leafs.end();
	while (beg != end) {
		MuCluster * x = *(beg++);
		vqueue.push(x);
	}
}
void _MSG_VSpace_down_top::visit_subsuming(MuCluster * x) {
	/* declarations */
	std::queue<MuCluster *> bfs_queue;
	std::set<MuCluster *> visit_set;

	/* initialization */
	bfs_queue.push(x); visit_set.insert(x);

	/* iterate by BFS from x to its parents or parents' parents */
	while (!bfs_queue.empty()) {
		/* get next node to be removed from VS */
		x = bfs_queue.front(); bfs_queue.pop();

		/* get edges from parents to x */
		const std::vector<MuSubsume> & edges
			= x->get_in_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while(beg != end) {
			/* get next unvisited parent */
			const MuSubsume & edge = *(beg++);
			MuCluster & parent = edge.get_source();
			if (visit_set.count(&parent) > 0) continue;

			/* set parent as visited */
			visitset.insert(&parent);

			/* update the bfs_queue */
			bfs_queue.push(&parent);
			visit_set.insert(&parent);
		}
	} /* end while: bfs_queue */

	/* return */ return;
}
bool _MSG_VSpace_down_top::accessible(MuCluster * x) {
	if (adset.count(x) == 0) return false;			/* not in the sub-graph */
	else if (visitset.count(x) > 0) return false;	/* have been visited */
	else {
		/* get the edges from x to its direct children */
		const std::vector<MuSubsume> & edges
			= x->get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();

		/* validate the x's children */
		while (beg != end) {
			/* get next child in sub-graph */
			const MuSubsume & edge = *(beg++);
			MuCluster & child = edge.get_target();
			if (adset.count(&child) == 0) continue;

			/* child is not visited yet */
			if (visitset.count(&child) == 0)
				return false;
		}
		return true;	/* all children have been visited */
	}
}
MuCluster * _MSG_VSpace_down_top::next() {
	MuCluster * next = nullptr;
	while (next == nullptr && !vqueue.empty()) {
		/* compute the next unvisited node in sub-graph */
		while (!vqueue.empty()) {
			next = vqueue.front(); 
			vqueue.pop(); vcache.push(next);

			if (visitset.count(next) == 0) {
				visitset.insert(next); 
				break;
			}
			else next = nullptr;
		} /* end while: vqueue */

		/* update vqueue from vcache */
		std::set<MuCluster *> records;
		while (!vcache.empty()) {
			MuCluster * x = vcache.front(); vcache.pop();

			const std::vector<MuSubsume> & edges
				= x->get_in_port().get_edges();
			auto beg = edges.begin(), end = edges.end();
			
			while (beg != end) {
				/* get next un-iterated parent */
				const MuSubsume & edge = *(beg++);
				MuCluster & parent = edge.get_source();
				if (records.count(&parent) > 0) continue;
				else records.insert(&parent);

				if (accessible(&parent))
					vqueue.push(&parent);
			}
		} /* end while: vcache */

	} /* end while */

	/* no more unvisited node */ return next;
}

_MSG_VSpace_top_down::~_MSG_VSpace_top_down() {
	while (!vqueue.empty()) vqueue.pop();
	while (!vcache.empty()) vcache.pop();
	visitset.clear();
}
void _MSG_VSpace_top_down::initial() {
	while (!vqueue.empty()) vqueue.pop();
	while (!vcache.empty()) vcache.pop();
	visitset.clear();

	auto beg = roots.begin();
	auto end = roots.end();
	while (beg != end) {
		MuCluster * x = *(beg++);
		vqueue.push(x);
	}
}
void _MSG_VSpace_top_down::visit_subsumed(MuCluster * x) {
	/* declarations */
	std::queue<MuCluster *> bfs_queue;
	std::set<MuCluster *> visit_set;

	/* initialization */
	bfs_queue.push(x); visit_set.insert(x);

	/* iterate by BFS from x to its parents or parents' parents */
	while (!bfs_queue.empty()) {
		/* get next node to be removed from VS */
		x = bfs_queue.front(); bfs_queue.pop();

		/* get edges from x to its children */
		const std::vector<MuSubsume> & edges
			= x->get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			/* get next unvisited parent */
			const MuSubsume & edge = *(beg++);
			MuCluster & child = edge.get_target();
			if (visit_set.count(&child) > 0) continue;

			/* set parent as visited */
			visitset.insert(&child);

			/* update the bfs_queue */
			bfs_queue.push(&child);
			visit_set.insert(&child);
		}
	} /* end while: bfs_queue */

	/* return */ return;
}
bool _MSG_VSpace_top_down::accessible(MuCluster * x) {
	if (adset.count(x) == 0) return false;		/* not in sub-graph */
	else if (visitset.count(x) > 0) return false;	/* have been visited */
	else {
		const std::vector<MuSubsume> & edges
			= x->get_in_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		
		while (beg != end) {
			/* get next parent in sub-graph */
			const MuSubsume & edge = *(beg++);
			MuCluster & parent = edge.get_source();
			if (adset.count(&parent) == 0) continue;

			/* not visited, then not accessible */
			if (visitset.count(&parent) == 0)
				return false;
		}
		return true;
	}
}
MuCluster * _MSG_VSpace_top_down::next() {
	MuCluster * next = nullptr;
	while (next == nullptr && !vqueue.empty()) {
		/* compute the next unvisited node in sub-graph */
		while (!vqueue.empty()) {
			next = vqueue.front();
			vqueue.pop(); vcache.push(next);

			if (visitset.count(next) == 0) {
				visitset.insert(next);
				break;
			}
			else next = nullptr;
		} /* end while: vqueue */

		  /* update vqueue from vcache */
		std::set<MuCluster *> records;
		while (!vcache.empty()) {
			MuCluster * x = vcache.front();
			vcache.pop();

			const std::vector<MuSubsume> & edges
				= x->get_ou_port().get_edges();
			auto beg = edges.begin(), end = edges.end();
			while (beg != end) {
				/* get next un-iterated parent */
				const MuSubsume & edge = *(beg++);
				MuCluster & child = edge.get_target();
				if (records.count(&child) > 0) continue;
				else records.insert(&child);

				if (accessible(&child))
					vqueue.push(&child);
			}
		} /* end while: vcache */

	} /* end while */

	/* no more unvisited node */ return next;
}

_MSG_VSpace_randomly::~_MSG_VSpace_randomly() {
	visitset.clear();
}
void _MSG_VSpace_randomly::initial() {
	visitset.clear();
	auto beg = adset.begin(), end = adset.end();
	while (beg != end) visitset.insert(*(beg++));
}
void _MSG_VSpace_randomly::visit_subsumed(MuCluster * x) {
	/* declarations */
	std::queue<MuCluster *> bfs_queue;
	std::set<MuCluster *> visit_set;

	/* initialization */
	bfs_queue.push(x); visit_set.insert(x);

	/* iterate by BFS from x to its parents or parents' parents */
	while (!bfs_queue.empty()) {
		/* get next node to be removed from VS */
		x = bfs_queue.front(); bfs_queue.pop();

		/* get edges from x to its children */
		const std::vector<MuSubsume> & edges
			= x->get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			/* get next unvisited parent */
			const MuSubsume & edge = *(beg++);
			MuCluster & child = edge.get_target();
			if (visit_set.count(&child) > 0) continue;

			/* set parent as visited */
			if(visitset.count(&child) > 0)
				visitset.erase(&child);

			/* update the bfs_queue */
			bfs_queue.push(&child);
			visit_set.insert(&child);
		}
	} /* end while: bfs_queue */

	/* return */ return;
}
void _MSG_VSpace_randomly::visit_subsuming(MuCluster * x) {
	/* declarations */
	std::queue<MuCluster *> bfs_queue;
	std::set<MuCluster *> visit_set;

	/* initialization */
	bfs_queue.push(x); visit_set.insert(x);

	/* iterate by BFS from x to its parents or parents' parents */
	while (!bfs_queue.empty()) {
		/* get next node to be removed from VS */
		x = bfs_queue.front(); bfs_queue.pop();

		/* get edges from parents to x */
		const std::vector<MuSubsume> & edges
			= x->get_in_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			/* get next unvisited parent */
			const MuSubsume & edge = *(beg++);
			MuCluster & parent = edge.get_source();
			if (visit_set.count(&parent) > 0) continue;

			/* set parent as visited */
			if (visitset.count(&parent) > 0)
				visitset.erase(&parent);

			/* update the bfs_queue */
			bfs_queue.push(&parent);
			visit_set.insert(&parent);
		}
	} /* end while: bfs_queue */

	/* return */ return;
}
MuCluster * _MSG_VSpace_randomly::next() {
	MuCluster * next = nullptr;
	if (!visitset.empty()) {
		auto beg = visitset.begin();
		next = *(beg);
		visitset.erase(next);
	}
	return next;
}

void MSGLinker::connect(MSGraph & g, OrderOption opt) {
	const MuHierarchy & hierarchy = g.get_hierarchy();
	int i, n = hierarchy.size_of_degress();
	std::map<MuCluster *, std::set<MuCluster *> *> solutions;

	this->open(g, opt);
	for (i = n - 1; i >= 0; i--) {
		/* get level at H[k] */
		const std::set<MuCluster *> & level 
			= hierarchy.get_clusters_at(i);
		auto beg = level.begin(), end = level.end();

		/* compute DS for each x in H[k] */
		while (beg != end) {
			MuCluster * x = *(beg++);
			std::set<MuCluster *> * DS = new std::set<MuCluster *>();
			compute_direct_subsumption(*x, *DS);
			solutions[x] = DS;
		}

		/* connect x to its DS */
		auto sbeg = solutions.begin();
		auto send = solutions.end();
		while (sbeg != send) {
			MuCluster * x = sbeg->first;
			std::set<MuCluster *> * DS = sbeg->second;
			connect_nodes(*x, *DS); delete DS; sbeg++;
		}
		solutions.clear();

		/* add nodes in H[k] to subgraph */
		add_nodes_in(level);
	}
	this->close();
}
void MSGLinker::open(MSGraph & g, OrderOption opt) {
	close(); graph = &g; g.clear_edges();
	switch (opt) {
	case down_top:
		vspace = new _MSG_VSpace_down_top(adset, leafs); break;
	case top_down:
		vspace = new _MSG_VSpace_top_down(adset, roots); break;
	case randomly:
		vspace = new _MSG_VSpace_randomly(adset); break;
	default:
		CError error(CErrorType::InvalidArguments, "MSGLinker::open", 
			"Unknown option (" + std::to_string(opt) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
}
void MSGLinker::close() {
	if (graph != nullptr) {
		graph->update_roots_and_leafs();
		adset.clear(); roots.clear(); leafs.clear();
		delete vspace; graph = nullptr;
	}
}
void MSGLinker::connect_nodes(MuCluster & x, 
	const std::set<MuCluster *> & DS) {
	if (DS.empty()) return;
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MuCluster & y = *(*(beg++));
		graph->connect(x, y);
	}
} 
void MSGLinker::add_nodes_in(const std::set<MuCluster *> & nodes) {
	std::set<MuCluster *> erase_roots;
	auto beg = roots.begin(), end = roots.end();
	while (beg != end) {
		MuCluster * x = *(beg++);
		if (x->get_in_port().get_degree() > 0)
			erase_roots.insert(x);
	}

	beg = erase_roots.begin(), end = erase_roots.end();
	while (beg != end) roots.erase(*(beg++));

	beg = nodes.begin(), end = nodes.end();
	while (beg != end) {
		MuCluster * x = *(beg++);
		adset.insert(x);

		if (x->get_ou_port().get_degree() == 0)
			leafs.insert(x);
		if (x->get_in_port().get_degree() == 0)
			roots.insert(x);
	}
}
void MSGLinker::compute_direct_subsumption(MuCluster & x, std::set<MuCluster *> & DS) {
	DS.clear(); vspace->initial(); MuCluster * y;

	while ((y = vspace->next()) != nullptr) {
		if (subsume(x, *y)) {
			DS.insert(y);
			vspace->visit_subsumed(y);
		}
		else vspace->visit_subsuming(y);
	}

	eliminate_DS(DS);	/* eliminate redundant mutants */
}
bool MSGLinker::subsume(MuCluster & x, MuCluster & y) {
	const BitSeq & xv = x.get_score_vector();
	const BitSeq & yv = y.get_score_vector();
	return xv.subsume(yv);
}
void MSGLinker::eliminate_DS(std::set<MuCluster *> & DS) {
	/* declarations */
	std::set<MuCluster *> eliminates;
	std::queue<MuCluster *> DS_queue;
	MuCluster * x, *y;

	/* initialize the queue for iteration */
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) DS_queue.push(*(beg++));

	/* eliminate all redundant mutants from DS */
	while (!DS_queue.empty()) {
		/* get next not-eliminated node */
		x = DS_queue.front(); DS_queue.pop();
		if (DS.count(x) == 0) continue;

		/* compute all those subsumed by x in current DS */
		beg = DS.begin(), end = DS.end();
		while (beg != end) {
			y = *(beg++);
			if (x == y) continue;
			else if (subsume(*x, *y))
				eliminates.insert(y);
		}

		/* eliminate them from DS */
		beg = eliminates.begin();
		end = eliminates.end();
		while (beg != end) DS.erase(*(beg++));
		eliminates.clear();
	}

	/* end */ return;
}