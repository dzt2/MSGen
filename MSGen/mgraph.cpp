#include "mgraph.h"
#include <algorithm>

unsigned int times;

MuCluster::MuCluster(MSGraph & g, MuCluster::ID id, const BitSeq & bits)
	: graph(g), cluster_id(id), score_vector(bits), score_degree(0) {
	MutantSpace & mspace = graph.get_space();
	mutants = mspace.create_set();
	in_port = new MuSubsumePort();
	ou_port = new MuSubsumePort();

	BitSeq::size_t i = 0, n = bits.bit_number();
	while (i < n) {
		if (score_vector.get_bit(i++) == BIT_1)
			score_degree++;
	}
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
	mutants->clear();

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
MuCluster * MSGraph::new_cluster(const BitSeq & bits) {
	MuCluster * cluster = new MuCluster(*this, clusters.size(), bits);
	clusters.push_back(cluster); hierarchy.add(*cluster);
	return cluster;
}
void MSGraph::add_mutant(MuCluster & cluster, Mutant::ID mid) {
	if (index.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments, 
			"MSGraph::add_mutant", 
			"Invalid mutant-id (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error); 
		exit(CErrorType::InvalidArguments);
	}
	else {
		mutants->add_mutant(mid);
		cluster.add_mutant(mid);
		index[mid] = &cluster;
	}
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

	/* efficiency analysis */ times = 0;

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

		times++;	/* efficiency analysis */
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

void MSGBuilder::open(MSGraph & g) {
	close(); graph = &g; g.clear();
	trie = new BitTrieTree();
}
void MSGBuilder::add(Mutant::ID mid, const BitSeq & bits) {
	if (graph == nullptr || trie == nullptr) {
		CError error(CErrorType::Runtime, "MSGBuilder::add", "Invalid access: not-opened");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}
	else {
		/* find the leaf where the cluster is referred by bit-string */
		BitTrie * leaf = trie->insert_vector(bits);

		/* the first time to create cluster and insert it to the trie */
		if (leaf->get_data() == nullptr) {
			MuCluster * cluster = graph->new_cluster(bits);
			graph->add_mutant(*cluster, mid); 
			leaf->set_data(cluster);
		}
		/* only to insert the mutant into the cluster */
		else {
			MuCluster * cluster = (MuCluster *)(leaf->get_data());
			graph->add_mutant(*cluster, mid);
		}
	}
}
void MSGBuilder::link() {
	graph->sort(); linker.connect(*graph);
}
void MSGBuilder::link(MSGLinker::OrderOption option) {
	graph->sort(); linker.connect(*graph, option);
}
void MSGBuilder::close() {
	if (graph != nullptr) {
		graph = nullptr;
		delete trie;
	}
}

void MSGraphPrinter::write(MSGraph & graph) {
	if (dir == nullptr) {
		CError error(CErrorType::Runtime, 
			"MSGraphPrinter::write(graph)", 
			"directory is not opened");
		CErrorConsumer::consume(error);
	}
	else {
		std::ofstream out1(dir->get_path() + "/graph.txt");
		write_mutant_graph(graph, out1); out1.close();

		std::ofstream out2(dir->get_path() + "/mutantLib.txt");
		write_mutant_lib(graph, out2); out2.close();
	}
}
void MSGraphPrinter::write_mutant_graph(MSGraph & graph, std::ostream & out) {
	/* Line : source |--> directly_subsumed_cluster(s) */
	size_t cnum = graph.size();
	out << "source\tdegree\tnext(s)\n";

	/* get each cluster and their edges */
	for (MuCluster::ID cid = 0; cid < cnum; cid++) {
		MuCluster & src = graph.get_cluster(cid);
		size_t degree = src.get_score_degree();

		out << cid << "\t" << degree << "\t";

		const std::vector<MuSubsume> & edges = src.get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			const MuSubsume & edge = *(beg++);
			out << edge.get_target().get_id();
			if (beg != end) out << "; ";
		}

		out << "\n";
	}

	/* return */ out << std::endl; return;
}
void MSGraphPrinter::write_mutant_lib(MSGraph & graph, std::ostream & out) {
	MutantSpace & mspace = graph.get_space();
	Mutant::ID mid = 0, mnum = mspace.number_of_mutants();

	out << "mid\toperator\torigin\treplace\tcluster\tdegree\n";

	for (mid = 0; mid < mnum; mid++) {
		if (graph.has_cluster_of(mid)) {
			Mutant & mutant = mspace.get_mutant(mid);
			MuCluster & cluster = graph.get_cluster_of(mid);
			
			out << mid << "\t" << mutant.get_operator() << "\t";

			const Mutation & mutation = 
				mutant.get_mutation(mutant.get_orders() - 1);
			std::string origin, replace;

			origin = mutation.get_location().get_text_at();
			replace = mutation.get_replacement();
			replace_text(origin); replace_text(replace);

			out << origin << "\t" << replace << "\t";

			out << cluster.get_id() << "\t" << cluster.get_score_degree() << "\n";
		}
	}
}
void MSGraphPrinter::replace_text(std::string & text) {
	int n = text.length();
	std::string replace = "";
	for (int i = 0; i < n; i++) {
		char ch = text[i];
		if (ch == '\t' || ch == '\n');
		else replace += ch;
	}
	text = replace;
}