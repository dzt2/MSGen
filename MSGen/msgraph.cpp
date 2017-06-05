#include "msgraph.h"
#include <algorithm>

MSGVertex::MSGVertex(const MSGraph & g) : graph(g), vid(-1), cluster(nullptr),
feature(nullptr), in_port(*this), ou_port(*this) {
	cluster = g.get_space().create_set();
}
MSGVertex::~MSGVertex() {
	graph.get_space().delete_set(cluster);
	if (feature != nullptr) delete feature;
}

MSGVertex & MSGraph::get_vertex(MSGVertex::ID vid) const {
	MSGVertex * vertex = nullptr;
	if (vid < vertices.size())
		vertex = vertices[vid];
	else {
		CError error(CErrorType::InvalidArguments, "MSGraph::get_vertex", "Invalid vid: " + std::to_string(vid));
		CErrorConsumer::consume(error);
	}
	return *vertex;
}
MSGVertex & MSGraph::new_vertex(const ScoreVector & svec) {
	MSGVertex * vex = new MSGVertex(*this);
	vex->set_feature(svec.get_vector(), svec.get_degree());
	vex->add_mutant(svec.get_mutant());
	vpool.insert(vex);
	return *vex;
}
bool MSGraph::add_cluster(MSGVertex & vex, const ScoreVector & svec) {
	vex.add_mutant(svec.get_mutant());
	return true;
}
bool MSGraph::connect(MSGVertex & x, MSGVertex & y) {
	x.link_to(y); y.link_by(x); return true;
}
bool MSGraph::add_vertex(MSGVertex & vex) {
	vex.set_id(vertices.size());
	vertices.push_back(&vex);
	return true;
}
void MSGraph::clear() {
	auto beg = vpool.begin(), end = vpool.end();
	while (beg != end) delete *(beg++);
	vpool.clear(); vertices.clear();
	roots.clear(); leafs.clear();
}

const std::vector<MSGVertex *> & MuHierarchy::get_vertices_at(BitSeq::size_t deg) const {
	const std::vector<MSGVertex *> * level = nullptr;
	if (deg_map.count(deg) > 0) {
		auto iter = deg_map.find(deg);
		level = iter->second;
	}
	else {
		CError error(CErrorType::InvalidArguments, "MuHierarchy::get_vertices_at", "Invalid deg: " + std::to_string(deg));
		CErrorConsumer::consume(error);
	}
	return *level;
}
bool MuHierarchy::add_vertex(MSGVertex & vex) {
	if (vex.get_feature() != nullptr) {
		/* get the degree of vertex */
		const MuFeature & feature = *(vex.get_feature());
		BitSeq::size_t degree = feature.get_degree();

		/* get the level of vertices in hierarchy */
		if (deg_map.count(degree) == 0) {
			deg_map[degree] = new std::vector<MSGVertex *>();
			deg_list.push_back(degree);
		}
		auto iter = deg_map.find(degree);
		auto level = iter->second;

		/* add the node in target level */
		level->push_back(&vex);
		return true;
	}
	else return false;
}
void MuHierarchy::sort_degrees() {
	std::sort(deg_list.begin(), deg_list.end());
}
void MuHierarchy::clear() {
	auto beg = deg_map.begin(), end = deg_map.end();
	while (beg != end) delete ((beg++)->second);
	deg_map.clear(); deg_list.clear();
}

void MSGIter_DownTop::open(const MSGraph & g) {
	close(); graph = &g;

	/* set the window at first */
	const std::set<MSGVertex *> & leafs = g.get_leafs();
	auto beg = leafs.begin(), end = leafs.end();
	while (beg != end)
		window.insert(*(beg++));
}
MSGVertex * MSGIter_DownTop::get_next() {
	if (!window.empty() || !qset.empty()) {
		/* load nodes from queue */
		if (window.empty()) {
			while (!qset.empty()) {
				/* get the next node in queue */
				MSGVertex & target = *(qset.front());
				qset.pop();

				/* get in-edges */
				const std::list<MuSubsume> & edges
					= target.get_in_port().get_edges();
				auto beg = edges.begin(), end = edges.end();
				while (beg != end) {
					const MSGVertex & source = (*(beg++)).get_source();
					MSGVertex * sptr = (MSGVertex *)(&source);
					if (accessible(sptr)) { window.insert(sptr); }
				}
			}
		}

		/* get the next for return */
		if (!window.empty()) {
			MSGVertex * vex = *(window.begin());
			window.erase(vex); visits.insert(vex);
			qset.push(vex); return vex;
		}
	}
	return nullptr;
}
void MSGIter_DownTop::erase(MSGVertex * vex) {
	if (visits.count(vex) == 0) {
		visits.insert(vex);
		/*if (window.count(vex) > 0) {
		window.erase(vex);
		qset.push(vex);
		}*/
	}
}
bool MSGIter_DownTop::accessible(MSGVertex * vex) {
	if (visits.count(vex) == 0) {
		if (graph->has_vertex(vex->get_id())) {
			/* get output edges */
			const std::list<MuSubsume> & edges
				= (vex->get_ou_port()).get_edges();
			auto beg = edges.begin(), end = edges.end();

			/* if any one of its child is not visited, it cannot be visited next */
			while (beg != end) {
				const MSGVertex & target = ((*(beg++)).get_target());
				if (!(graph->has_vertex(target.get_id()))) continue;
				MSGVertex * tptr = (MSGVertex *)(&target);
				if (visits.count(tptr) == 0) return false;
			}
			return true;
		}
	}
	return false;
}

void MSGIter_TopDown::open(const MSGraph & g) {
	close(); graph = &g;

	/* set the window at first */
	const std::set<MSGVertex *> & roots = g.get_roots();
	auto beg = roots.begin(), end = roots.end();
	while (beg != end)
		window.insert(*(beg++));
}
MSGVertex * MSGIter_TopDown::get_next() {
	if (!window.empty() || !qset.empty()) {
		/* load nodes from queue */
		if (window.empty()) {
			while (!qset.empty()) {
				/* get the next node in queue */
				MSGVertex & source = *(qset.front());
				qset.pop();

				/* get in-edges */
				const std::list<MuSubsume> & edges
					= source.get_ou_port().get_edges();
				auto beg = edges.begin(), end = edges.end();
				while (beg != end) {
					const MSGVertex & target = (*(beg++)).get_target();
					MSGVertex * tptr = (MSGVertex *)(&target);
					if (accessible(tptr)) { window.insert(tptr); }
				}
			}
		}

		/* get the next for return */
		if (!window.empty()) {
			MSGVertex * vex = *(window.begin());
			window.erase(vex); qset.push(vex);
			visits.insert(vex); return vex;
		}
	}
	return nullptr;
}
void MSGIter_TopDown::erase(MSGVertex * vex) {
	if (visits.count(vex) == 0) {
		visits.insert(vex);
		/*if (window.count(vex) > 0) {
		window.erase(vex);
		qset.push(vex);
		}*/
	}
}
bool MSGIter_TopDown::accessible(MSGVertex * vex) {
	if (visits.count(vex) == 0) {
		if (graph->has_vertex(vex->get_id())) {
			/* get input edges */
			const std::list<MuSubsume> & edges
				= (vex->get_in_port()).get_edges();
			auto beg = edges.begin(), end = edges.end();

			/* if any one of its parent is not visited, it cannot be visited next */
			while (beg != end) {
				const MSGVertex & source = ((*(beg++)).get_source());
				if (!(graph->has_vertex(source.get_id()))) continue;
				MSGVertex * sptr = (MSGVertex *)(&source);
				if (visits.count(sptr) == 0) return false;
			}
			return true;
		}
	}
	return false;
}

void MSGIter_Random::open(const MSGraph & g) {
	close(); graph = &g;

	/* set random sequence */
	size_t k = 0, n = g.number_of_vertices();
	while (k < n) {
		MSGVertex & vex = g.get_vertex(k++);
		nodes.insert(&vex);
	}
}
MSGVertex * MSGIter_Random::get_next() {
	if (!nodes.empty()) {
		MSGVertex * vex = *(nodes.begin());
		nodes.erase(vex); return vex;
	}
	else return nullptr;
}
void MSGIter_Random::erase(MSGVertex * vex) {
	nodes.erase(vex);
}

bool MSGBuilder::open(MSGraph & g) {
	close();	/* remove original state */

				/* clear graph */
	graph = &g; graph->clear();
	/* clear hierarchy */
	hierarchy.clear();
	/* initialize trie */
	trie = new BitTrieTree();

	return true;
}
bool MSGBuilder::clusterMutants(ScoreProducer & producer, ScoreConsumer & consumer) {
	ScoreVector * svector;
	while ((svector = producer.produce()) != nullptr) {
		/* insert the vector into bit-trie */
		BitTrie * leaf = trie->insert_vector(svector->get_vector());

		/* first created, then construct a vertex */
		if (leaf->get_data() == nullptr) {
			MSGVertex & vex = graph->new_vertex(*svector);
			hierarchy.add_vertex(vex);
			leaf->set_data(&vex);
		}
		/* existing node, then add mutant */
		else {
			MSGVertex * vex = (MSGVertex *)(leaf->get_data());
			graph->add_cluster(*vex, *svector);
		}

		/* delete allocated score vector */
		consumer.consume(svector);
	}
	return true;
}
bool MSGBuilder::rankMutantByDegree() {
	hierarchy.sort_degrees();
	return true;
}
bool MSGBuilder::build_graph() {
	/* get the hierarchy levels */
	const std::vector<BitSeq::size_t> &
		deg_list = hierarchy.get_degrees();
	int k = deg_list.size() - 1;
	std::set<MSGVertex *> DS; size_t E = 0;

	/* build the graph from leafs to roots */
	while (k >= 0) {
		/* get H[k] */
		BitSeq::size_t degree = deg_list[k];
		const std::vector<MSGVertex *> & level
			= hierarchy.get_vertices_at(degree);

		/* link x to the graph */
		auto lbeg = level.begin(), lend = level.end(); int c = 0;
		while (lbeg != lend) {
			MSGVertex & x = *(*(lbeg++));	/* get x in H[k] */
			direct_subsumed(x, DS);	/* calculate DS */
			E += DS.size(); c++;

			auto dbeg = DS.begin(), dend = DS.end();
			while (dbeg != dend) {
				MSGVertex & y = *(*(dbeg++));
				graph->connect(x, y);
			}
		}

		/*add H[k] into graph*/ add_vertices(level);

		/* to H[k-1] continues */ k = k - 1;
	}

	return true;
}
void MSGBuilder::direct_subsumed(const MSGVertex & x, std::set<MSGVertex *> & DS) {
	/* choose a iterator for implement */
	//direct_subsumed_downtop(x, DS);
	direct_subsumed_topdown(x, DS);
	//direct_subsumed_random(x, DS);
}
void MSGBuilder::add_vertices(const std::vector<MSGVertex *> & level) {
	/* update original roots, leafs are not needed to be updated */
	std::set<MSGVertex *> not_roots;
	auto rbeg = (graph->roots).begin();
	auto rend = (graph->roots).end();
	while (rbeg != rend) {
		MSGVertex & node = *(*(rbeg++));
		if (node.get_in_degree() > 0)
			not_roots.insert(&node);
	}

	/* remove invalid roots from graph */
	auto rm_beg = not_roots.begin();
	auto rm_end = not_roots.end();
	while (rm_beg != rm_end) {
		(graph->roots).erase(*(rm_beg++));
	}

	/* add both roots and leafs */
	auto beg = level.begin(), end = level.end();
	while (beg != end) {
		MSGVertex & vex = *(*(beg++));
		graph->add_vertex(vex);
		if (vex.get_in_degree() == 0)
			(graph->roots).insert(&vex);
		if (vex.get_ou_degree() == 0)
			(graph->leafs).insert(&vex);
	}
}

void MSGBuilder::direct_subsumed_downtop(const MSGVertex & x, std::set<MSGVertex *> & DS) {
	/* initialization */
	DS.clear(); dt_iter.open(x.get_graph());
	MSGVertex * y; std::set<MSGVertex *> rm;

	/* iterate to compute DS */
	while ((y = dt_iter.get_next()) != nullptr) {
		if (is_subsume(x, *y)) {
			subsumed(*y, rm);
			auto beg = rm.begin(), end = rm.end();
			while (beg != end) {
				DS.erase(*beg);
				//dt_iter.erase(*beg);
				beg++;
			}
			DS.insert(y);
		}
		else {
			subsuming(*y, rm);
			auto beg = rm.begin(), end = rm.end();
			while (beg != end) {
				dt_iter.erase(*beg); beg++;
			}
		}
	}
}
void MSGBuilder::direct_subsumed_topdown(const MSGVertex & x, std::set<MSGVertex *> & DS) {
	/* initialization */
	DS.clear(); td_iter.open(x.get_graph());
	MSGVertex * y; std::set<MSGVertex *> rm;

	/* iterate node from visit space */
	while ((y = td_iter.get_next()) != nullptr) {
		if (is_subsume(x, *y)) {
			subsumed(*y, rm);
			auto beg = rm.begin(), end = rm.end();
			while (beg != end) {
				DS.erase(*beg);
				td_iter.erase(*beg);
				beg++;
			}
			DS.insert(y);
		}
		/*else {
		subsuming(*y, rm);
		auto beg = rm.begin(), end = rm.end();
		while (beg != end) {
		td_iter.erase(*beg); beg++;
		}
		}*/
	}
}
void MSGBuilder::direct_subsumed_random(const MSGVertex & x, std::set<MSGVertex *> & DS) {
	/* initialization */
	DS.clear(); rd_iter.open(x.get_graph());
	MSGVertex * y; std::set<MSGVertex *> rm;

	/* iterate node from visit space */
	while ((y = rd_iter.get_next()) != nullptr) {
		if (is_subsume(x, *y)) {
			subsumed(*y, rm);
			auto beg = rm.begin(), end = rm.end();
			while (beg != end) {
				DS.erase(*beg);
				rd_iter.erase(*beg);
				beg++;
			}
			DS.insert(y);
		}
		else {
			subsuming(*y, rm);
			auto beg = rm.begin(), end = rm.end();
			while (beg != end) {
				rd_iter.erase(*beg); beg++;
			}
		}
	}
}

bool MSGBuilder::subsumed(const MSGVertex & x, std::set<MSGVertex *> & set) {
	/* initialization */
	std::queue<MSGVertex *> qlist;
	qlist.push((MSGVertex *)(&x));

	set.clear();
	while (!qlist.empty()) {
		/* get next nodes in queue */
		MSGVertex * z = qlist.front(); qlist.pop();

		/* get ouput edges */
		const MSGVexPort & port = z->get_ou_port();
		const std::list<MuSubsume> &edges = port.get_edges();

		/* iterate directly subsumed ones in graph */
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			const MSGVertex & target = (*(beg++)).get_target();
			MSGVertex * tptr = (MSGVertex *)(&target);

			if (set.count(tptr) == 0 && graph->has_vertex(target.get_id())) {
				set.insert(tptr); qlist.push(tptr);
			}
		}
	} /* end while: qlist */

	/* return */ return true;
}
bool MSGBuilder::subsuming(const MSGVertex & x, std::set<MSGVertex *> & set) {
	/* initialization */
	std::queue<MSGVertex *> qlist;
	qlist.push((MSGVertex *)(&x));

	set.clear();
	while (!qlist.empty()) {
		/* get next nodes in queue */
		MSGVertex * z = qlist.front(); qlist.pop();

		/* get ouput edges */
		const MSGVexPort & port = z->get_in_port();
		const std::list<MuSubsume> &edges = port.get_edges();

		/* iterate directly subsumed ones in graph */
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			const MSGVertex & source = (*(beg++)).get_source();
			MSGVertex * sptr = (MSGVertex *)(&source);

			if (set.count(sptr) == 0 && graph->has_vertex(source.get_id())) {
				set.insert(sptr); qlist.push(sptr);
			}
		}
	} /* end while: qlist */

	/* return */ return true;
}
bool MSGBuilder::is_subsume(const MSGVertex & x, const MSGVertex & y) {
	const BitSeq & xvec = x.get_feature()->get_vector();
	const BitSeq & yvec = y.get_feature()->get_vector();
	return xvec.subsume(yvec);
}

bool MSGBuilder::build(MSGraph & graph, ScoreProducer & producer, ScoreConsumer & consumer) {
	this->open(graph);
	this->clusterMutants(producer, consumer);
	this->rankMutantByDegree();
	this->build_graph();
	return true;
}
void MSGBuilder::add_score_vector(const ScoreVector & svec) {
	/* insert the vector into bit-trie */
	BitTrie * leaf = trie->insert_vector(svec.get_vector());

	/* first created, then construct a vertex */
	if (leaf->get_data() == nullptr) {
		MSGVertex & vex = graph->new_vertex(svec);
		hierarchy.add_vertex(vex);
		leaf->set_data(&vex);
	}
	/* existing node, then add mutant */
	else {
		MSGVertex * vex = (MSGVertex *)(leaf->get_data());
		graph->add_cluster(*vex, svec);
	}
}
void MSGBuilder::end_score_vectors() {
	this->rankMutantByDegree();
	this->build_graph();
}