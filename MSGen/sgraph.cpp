#include "sgraph.h"
#include <algorithm>

MSG_Port::~MSG_Port() {
	int n = edges.size(), k;
	for (k = 0; k < n; k++) 
		delete edges[k];
	edges.clear();
}
MSG_Edge & MSG_Port::get_edge(int k) const {
	if (k < 0 || k >= edges.size()) {
		CError error(CErrorType::InvalidArguments,"MSG_Port::get_edge(k)","Invalid k: " + std::to_string(k));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *edges[k];
}
bool MSG_Port::link(MSG_Node & x, MSG_Node & y) {
	MSG_Edge * edge = new MSG_Edge(x, y);
	edges.push_back(edge); return true;
}

MSG_Node::MSG_Node(MS_Graph & g, long cid, const BitSeq & svec) : 
	graph(g), id(cid), mutants(nullptr), score_vector(svec), in_port(), ou_port() {
	mutants = g.get_space().create_set();
	score_degree = svec.degree();
}
MSG_Node::~MSG_Node() {
	graph.get_space().delete_set(mutants);
}
bool MSG_Node::link_to(MSG_Node & next) {
	ou_port.link(*this, next);
	next.in_port.link(*this, next);
	return true;
}

MSG_Node & MS_Graph::get_node(long id) const {
	if (id < 0 || id >= nodes.size()) {
		CError error(CErrorType::InvalidArguments, "MS_Graph::get_node(k)", "Invalid k: " + std::to_string(id));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *nodes[id];
}
MSG_Node & MS_Graph::get_node_of(Mutant::ID mid) const {
	if (mut_node.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, "MS_Graph::get_node_of(mid)", 
			"Undefined mid: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = mut_node.find(mid); return *(iter->second);
	}
}
MSG_Node & MS_Graph::new_node(const BitSeq & svec) {
	MSG_Node * node = new MSG_Node(*this, nodes.size(), svec);
	nodes.push_back(node); return *node;
}
bool MS_Graph::add_mutant(MSG_Node & node, Mutant::ID mid) {
	if (mut_node.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments, "MS_Graph::add_mutant(node, mid)", 
			"Invalid mid: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		node.get_mutants().add_mutant(mid);
		mut_node[mid] = &node; return true;
	}
}
bool MS_Graph::clear() {
	int k, n = nodes.size();
	for (k = 0; k < n; k++) 
		delete nodes[k];
	nodes.clear(); mut_node.clear();
	return true;
}

void MSG_Build::open(ScoreProducer & producer, ScoreConsumer & consumer) {
	close(); graph.clear();
	this->producer = &producer;
	this->consumer = &consumer;
}
void MSG_Build::close() {
	producer = nullptr;
	consumer = nullptr;
}
bool MSG_Build::build() {
	if (producer == nullptr || consumer == nullptr) {
		CError error(CErrorType::InvalidArguments, "MSG_Build::build()", "not openned yet");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return construct();
}

void MSG_Build_Exhaustive::load_vectors() {
	this->clear_vectors();
	ScoreVector * score_vector;
	while ((score_vector = producer->produce()) != nullptr)
		vectors.push_back(score_vector);

	size = vectors.size();
	matrix = new bool *[size];
	for (int i = 0; i < size; i++) {
		matrix[i] = new bool[size];
		for (int j = 0; j < size; j++)
			matrix[i][j] = false;
	}

	return;
}
void MSG_Build_Exhaustive::clear_vectors() {
	if (vectors.empty()) return;

	int i, n = vectors.size();
	for (i = 0; i < n; i++)
		consumer->consume(vectors[i]);
	vectors.clear(); 
	
	for (int i = 0; i < size; i++)
		delete matrix[i];
	delete matrix; size = 0;
}
bool MSG_Build_Exhaustive::construct() {
	this->load_vectors();

	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if (i == j) {
				matrix[i][j] = true;
			}
			else {
				ScoreVector & x = *(vectors[i]);
				ScoreVector & y = *(vectors[j]);
				if (this->subsume(x.get_vector(), y.get_vector()))
					matrix[i][j] = true;
			}
		}
	}

	this->clear_vectors();
	return true;
}

bool MSG_Build_Classical::classify_vectors(const std::set<ScoreVector *> & source,
	std::set<ScoreVector *> & indistinguishables, std::set<ScoreVector *> & subsumes,
	std::set<ScoreVector *> & subsumeds, std::set<ScoreVector *> & non_subsumes) {
	/* initialize the outputs */
	indistinguishables.clear(); non_subsumes.clear();
	subsumes.clear(); subsumeds.clear();
	if (source.empty()) return false;

	/* indistinguishable construct */
	auto beg = source.begin();
	auto end = source.end();
	ScoreVector & x = *(*(beg++));
	indistinguishables.insert(&x);

	/* classification */
	while (beg != end) {
		ScoreVector & y = *(*(beg++));
		if (subsume(x.get_vector(), y.get_vector())) {
			if (subsume(y.get_vector(), x.get_vector())) 
				indistinguishables.insert(&y);
			else subsumeds.insert(&y);
		}
		else if (subsume(y.get_vector(), x.get_vector()))
			subsumes.insert(&y);
		else non_subsumes.insert(&y);
	}

	/* return */ return true;
}
MSG_Node & MSG_Build_Classical::build_up_node(const std::set<ScoreVector *> & indistinguishables) {
	/* create the node for indistinguishable mutants */
	ScoreVector * next = *(indistinguishables.begin());
	const BitSeq & bitseq = next->get_vector();
	MSG_Node & node = graph.new_node(bitseq);

	/* add these mutants to the node */
	auto beg = indistinguishables.begin();
	auto end = indistinguishables.end();
	while (beg != end) {
		ScoreVector & score_vector = *(*(beg++));
		graph.add_mutant(node, score_vector.get_mutant());
	}

	/* return */ return node;
}
bool MSG_Build_Classical::derive_roots_leafs(const std::set<ScoreVector *> & source,
	std::set<MSG_Node *> & roots, std::set<MSG_Node *> & leafs) {
	roots.clear(); leafs.clear();

	auto beg = source.begin();
	auto end = source.end();
	while (beg != end) {
		ScoreVector & next = *(*(beg++));
		MSG_Node & node = graph.get_node_of(next.get_mutant());
		if (node.get_in_port().degree() == 0) roots.insert(&node);
		if (node.get_ou_port().degree() == 0) leafs.insert(&node);
	}

	return true;
}
bool MSG_Build_Classical::topdown_ready(MSG_Node & next, const std::set<MSG_Node *> & visits) {
	const MSG_Port & port = next.get_in_port();
	for (int i = 0; i < port.degree(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & prev = edge.get_source();
		if (visits.count(&prev) == 0) return false;
	}
	return true;
}
bool MSG_Build_Classical::derive_direct_subsumed(MSG_Node & source, 
	const std::set<MSG_Node *> &troots, std::set<MSG_Node *> & DS) {
	/* initialization */
	std::queue<MSG_Node *> tqueue;
	std::set<MSG_Node *> tvisits, records;
	auto tbeg = troots.begin();
	auto tend = troots.end();
	while (tbeg != tend) {
		records.insert(*tbeg);
		tqueue.push(*(tbeg++));
	}

	DS.clear(); 
	while (!tqueue.empty()) {
		/* get the next node in queue */
		MSG_Node * target = tqueue.front();
		tqueue.pop(); tvisits.insert(target);

		if (subsume(source, *target)) 
			DS.insert(target);
		/* put the target's "valid" child to queue 
				1) not in queue yet;
				2) all the parents have been visited.
		*/
		else {
			const MSG_Port & port = target->get_ou_port();
			for (int i = 0; i < port.degree(); i++) {
				MSG_Edge & edge = port.get_edge(i);
				MSG_Node & next = edge.get_target();
				if (records.count(&next) == 0 &&
					topdown_ready(next, tvisits)) {
					tqueue.push(&next);
					records.insert(&next);
				}
			}
		}
	}

	return true;
}
bool MSG_Build_Classical::update_direct_subsumed(MSG_Node & source, std::set<MSG_Node *> & DS) {
	std::set<MSG_Node *> IS;
	const MSG_Port & port = source.get_ou_port();
	for (int i = 0; i < port.degree() && !DS.empty(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & next = edge.get_target();

		auto beg = DS.begin();
		auto end = DS.end();
		while (beg != end) {
			MSG_Node & check = *(*(beg++));
			if (subsume(next, check)) {
				IS.insert(&check);
			}
		}

		beg = IS.begin(), end = IS.end();
		while (beg != end) DS.erase(*(beg++));
		IS.clear();
	}
	return purify_direct_subsumed(DS);
}
bool MSG_Build_Classical::purify_direct_subsumed(std::set<MSG_Node *> & DS) {
	std::set<MSG_Node *> trash;
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MSG_Node * src = *(beg++);
		if (trash.count(src) > 0) continue;
		else {
			auto beg2 = DS.begin();
			while (beg2 != end) {
				MSG_Node * trg = *(beg2++);
				if (src == trg) continue;
				else if (trash.count(trg) > 0) continue;
				else if (subsume(*src, *trg))
					trash.insert(trg);
			}
		}
	}

	beg = trash.begin(), end = trash.end();
	while (beg != end) DS.erase(*(beg++));
	return true;
}
bool MSG_Build_Classical::combine(const std::set<MSG_Node *> & ALeafs,
	const std::set<MSG_Node *> & BRoots, std::map<MSG_Node *, std::set<MSG_Node *> *> & ans) {
	/* initialization */
	std::queue<MSG_Node *> squeue;
	std::set<MSG_Node *> svisits;
	auto sbeg = ALeafs.begin();
	auto send = ALeafs.end();
	while (sbeg != send) {
		squeue.push(*(sbeg));
		svisits.insert(*sbeg);
		sbeg++;
	}

	/* iterate */
	ans.clear(); std::set<MSG_Node *> DS;
	while (!squeue.empty()) {
		/* get next node in A-graph */
		MSG_Node & source = *(squeue.front()); squeue.pop();

		/* compute its direct subsumed nodes */
		derive_direct_subsumed(source, BRoots, DS);
		update_direct_subsumed(source, DS);

		/* record the answer for source */
		if (!DS.empty()) 
			ans[&source] = new std::set<MSG_Node *>(DS);

		/* to the upper level */
		const MSG_Port & port = source.get_in_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & prev = edge.get_source();
			if (svisits.count(&prev) == 0) {
				svisits.insert(&prev);
				squeue.push(&prev);
			}
		}
	}

	DS.clear(); return true;
}
bool MSG_Build_Classical::build_subgraph(const std::set<ScoreVector *> & origin) {
	if (origin.empty()) return false;

	/* declarations */
	std::set<ScoreVector *> I, S, D, N;
	std::set<MSG_Node *> SRoots, SLeafs;
	std::set<MSG_Node *> DRoots, DLeafs;
	std::set<MSG_Node *> NRoots, NLeafs;

	/* classify the first node and groups */
	classify_vectors(origin, I, S, D, N);
	/* construct node for indistinguishable ones */
	MSG_Node & node = build_up_node(I);

	/* compute graphs for S, D, N */
	build_subgraph(S);
	build_subgraph(D);
	build_subgraph(N);

	/* derive roots and leafs */
	derive_roots_leafs(S, SRoots, SLeafs);
	derive_roots_leafs(D, DRoots, DLeafs);
	derive_roots_leafs(N, NRoots, NLeafs);

	/* combine-subgraphs */
	std::map<MSG_Node *, std::set<MSG_Node *> *> ND, SN;
	combine(NLeafs, DRoots, ND); 
	combine(SLeafs, NRoots, SN);

	/* combine the node with sub-graph */
	auto sbeg = SLeafs.begin();
	auto send = SLeafs.end();
	while (sbeg != send) {
		MSG_Node & src = *(*(sbeg++));
		graph.connect(src, node);
	}
	auto dbeg = DRoots.begin();
	auto dend = DRoots.end();
	while (dbeg != dend) {
		MSG_Node & trg = *(*(dbeg++));
		graph.connect(node, trg);
	}

	/* combine N-->D */
	auto ND_beg = ND.begin();
	auto ND_end = ND.end();
	while (ND_beg != ND_end) {
		MSG_Node & src = *((*ND_beg).first);
		std::set<MSG_Node *> * nexts = (*ND_beg).second;

		auto next_beg = nexts->begin();
		auto next_end = nexts->end();
		while (next_beg != next_end) 
			graph.connect(src, *(*(next_beg++)));
		
		delete nexts; ND_beg++;
	}
	/* combine S-->N */
	auto SN_beg = SN.begin();
	auto SN_end = SN.end();
	while (SN_beg != SN_end) {
		MSG_Node & src = *((*SN_beg).first);
		std::set<MSG_Node *> * nexts = (*SN_beg).second;

		auto next_beg = nexts->begin();
		auto next_end = nexts->end();
		while (next_beg != next_end)
			graph.connect(src, *(*(next_beg++)));

		delete nexts; SN_beg++;
	}

	/* return */ return true;
}
bool MSG_Build_Classical::construct() {
	/* derive the score vector into building */
	std::set<ScoreVector *> vectors;
	ScoreVector * score_vector;
	while ((score_vector = producer->produce()) != nullptr) {
		vectors.insert(score_vector);
	}

	/* build up MSG */ build_subgraph(vectors);

	/* consume all the vectors */
	auto beg = vectors.begin();
	auto end = vectors.end();
	while (beg != end) 
		consumer->consume(*(beg++));

	/* return */ return true;
}

bool MSG_Build_Fast::clustering() {
	BitTrieTree trie; ScoreVector * vec;

	while ((vec = producer->produce()) != nullptr) {
		Mutant::ID mid = vec->get_mutant();
		const BitSeq & bits = vec->get_vector();

		BitTrie * leaf = trie.insert_vector(bits);
		if (leaf->get_data() == nullptr) {
			MSG_Node & node = graph.new_node(bits);
			leaf->set_data(&node); 
			graph.add_mutant(node, mid);
		}
		else {
			MSG_Node & node = *((MSG_Node *)(leaf->get_data()));
			graph.add_mutant(node, mid);
		}

		consumer->consume(vec);
	}

	return true;
}
bool MSG_Build_Fast::rankByDegree(std::vector<std::set<MSG_Node *> *> & H) {
	std::map<unsigned int, std::set<MSG_Node *> *> hmap;
	std::vector<unsigned> keys;
	for (int i = 0; i < graph.size(); i++) {
		MSG_Node & node = graph.get_node(i);
		unsigned degree = node.get_score_degree();
		if (hmap.count(degree) == 0) {
			hmap[degree] = new std::set<MSG_Node *>();
			keys.push_back(degree);
		}

		auto iter = hmap.find(degree);
		std::set<MSG_Node *> * level = iter->second;
		level->insert(&node);
	}

	std::sort(keys.begin(), keys.end());

	H.clear();
	for (int i = 0; i < keys.size(); i++) {
		unsigned degree = keys[i];
		auto iter = hmap.find(degree);
		H.push_back(iter->second);
	}
	
	return true;
}
bool MSG_Build_Fast::directSubsumed(MSG_Node & x, const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & DS) {
	return directSubsumed_topdown(x, VS, DS);
	//return directSubsumed_downtop(x, VS, DS);
	//return directSubsumed_randomy(x, VS, DS);
}
bool MSG_Build_Fast::construct() {
	// step1. clustering 
	clustering();

	// step2. ranking 
	std::vector<std::set<MSG_Node *> *> H;
	rankByDegree(H);

	// step3. linking from H[n] to H[1]
	std::set<MSG_Node *> base, * DS;
	std::map<MSG_Node *, std::set<MSG_Node *> *> ans;
	for (int k = H.size() - 1; k >= 0; k--) {
		/* get H[k] */
		std::set<MSG_Node *> & level = *(H[k]);

		/* record DS for each node in H[k] */
		auto beg = level.begin(), end = level.end();
		while (beg != end) {
			/* get the nodes in DS */
			MSG_Node & x = *(*(beg++));
			DS = new std::set<MSG_Node *>();
			directSubsumed(x, base, *DS); 
			ans[&x] = DS;
		}

		/* update the nodes and edges "in" DMSG */
		beg = level.begin(), end = level.end();
		while (beg != end) {
			/* add nodes in the DMSG */
			MSG_Node * x = *(beg++);
			base.insert(x);

			/* connect x to its DS nodes */
			auto iter = ans.find(x);
			DS = iter->second;
			auto dbeg = DS->begin();
			auto dend = DS->end();
			while (dbeg != dend) {
				MSG_Node & y = *(*(dbeg++));
				graph.connect(*x, y);
			}
			delete DS; DS = nullptr;
		}

		/* delete current hierarchy */ 
		delete &level; ans.clear();
	}

	// step4. end and delete resources.
	base.clear(); H.clear(); return true;
}

bool MSG_Build_Fast::derive_roots(const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & roots) {
	roots.clear();
	auto beg = VS.begin();
	auto end = VS.end();
	while (beg != end) {
		MSG_Node & x = *(*(beg++));
		if (x.get_in_port().degree() == 0) 
			roots.insert(&x);
	}
	return true;
}
bool MSG_Build_Fast::derive_leafs(const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & leafs) {
	leafs.clear();
	auto beg = VS.begin();
	auto end = VS.end();
	while (beg != end) {
		MSG_Node & x = *(*(beg++));
		if (x.get_ou_port().degree() == 0)
			leafs.insert(&x);
	}
	return true;
}
bool MSG_Build_Fast::purify_subsumeds(std::set<MSG_Node *> & DS) {
	std::set<MSG_Node *> trash;
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MSG_Node * src = *(beg++);
		if (trash.count(src) > 0) continue;
		else {
			auto beg2 = DS.begin();
			while (beg2 != end) {
				MSG_Node * trg = *(beg2++);
				if (src == trg) continue;
				else if (trash.count(trg) > 0) continue;
				else if (subsume(*src, *trg))
					trash.insert(trg);
			}
		}
	}

	beg = trash.begin(), end = trash.end();
	while (beg != end) DS.erase(*(beg++));
	return true;
}

bool MSG_Build_Fast::available_topdown(MSG_Node & y, const std::set<MSG_Node *> & VS) {
	const MSG_Port & port = y.get_in_port();
	for (int i = 0; i < port.degree(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & prev = edge.get_source();
		if (VS.count(&prev) > 0) return false;
	}
	return true;
}
bool MSG_Build_Fast::directSubsumed_topdown(MSG_Node & x, const std::set<MSG_Node *> & base, std::set<MSG_Node *> & DS) {
	/* get the roots in current DMSG set VS */
	std::set<MSG_Node *> roots; 
	derive_roots(base, roots);

	/* initialize the sequence to access nodes in VS */
	std::queue<MSG_Node *> queue;
	std::set<MSG_Node *> records;
	auto beg = roots.begin(), end = roots.end();
	while (beg != end) {
		queue.push(*(beg));
		records.insert(*(beg++));
	}

	/* iterate */
	DS.clear();
	std::set<MSG_Node *> VS(base);
	while (!queue.empty()) {
		/* get next accessible node */
		MSG_Node & y = *(queue.front()); 
		queue.pop();

		/* update the visit space */
		if (VS.count(&y) == 0) 
			continue;
		else VS.erase(&y); 
		
		/* x > y */
		if (subsume(x, y)) {
			DS.insert(&y);
			/* erase children for one-level (performance) */
			const MSG_Port & port = y.get_ou_port();
			for (int i = 0; i < port.degree(); i++) {
				MSG_Edge & edge = port.get_edge(i);
				MSG_Node & next = edge.get_target();
				VS.erase(&next);
			}
		}
		else {
			/* push children */
			const MSG_Port & port = y.get_ou_port();
			for (int i = 0; i < port.degree(); i++) {
				MSG_Edge & edge = port.get_edge(i);
				MSG_Node & next = edge.get_target();
				if (records.count(&next) == 0
					&& VS.count(&next) > 0
					&& available_topdown(next, VS)) {
					records.insert(&next);
					queue.push(&next);
				}
			}
		}
	}

	/* purify DS */ return purify_subsumeds(DS);
}
bool MSG_Build_Fast::available_downtop(MSG_Node & y, const std::set<MSG_Node *> & VS) {
	const MSG_Port & port = y.get_ou_port();
	for (int i = 0; i < port.degree(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & next = edge.get_target();
		if (VS.count(&next) > 0) return false;
	}
	return true;
}
bool MSG_Build_Fast::directSubsumed_downtop(MSG_Node & x, const std::set<MSG_Node *> & base, std::set<MSG_Node *> & DS) {
	/* get the leafs in current DMSG set VS */
	std::set<MSG_Node *> leafs;
	derive_leafs(base, leafs);

	/* initialize the sequence to access nodes in VS */
	std::queue<MSG_Node *> queue;
	std::set<MSG_Node *> records;
	auto beg = leafs.begin(), end = leafs.end();
	while (beg != end) {
		queue.push(*(beg));
		records.insert(*(beg++));
	}

	/* iterate */
	DS.clear();
	std::set<MSG_Node *> VS(base);
	while (!queue.empty()) {
		/* get next accessible node */
		MSG_Node & y = *(queue.front());
		queue.pop();

		/* update the visit space */
		if (VS.count(&y) == 0) continue;
		else VS.erase(&y);

		/* x > y */
		if (subsume(x, y)) {
			DS.insert(&y);
			/* push children */
			const MSG_Port & port = y.get_in_port();
			for (int i = 0; i < port.degree(); i++) {
				MSG_Edge & edge = port.get_edge(i);
				MSG_Node & prev = edge.get_source();
				if (records.count(&prev) == 0
					&& VS.count(&prev) > 0
					&& available_downtop(prev, VS)) {
					records.insert(&prev);
					queue.push(&prev);
				}
			}
		}
		else {
			const MSG_Port & port = y.get_in_port();
			for (int i = 0; i < port.degree(); i++) {
				MSG_Edge & edge = port.get_edge(i);
				MSG_Node & prev = edge.get_source();
				VS.erase(&prev);
			}
		}
	}

	/* purify DS */ return purify_subsumeds(DS);
}
bool MSG_Build_Fast::erase_subsuming(MSG_Node & y, std::set<MSG_Node *> & VS) {
	std::queue<MSG_Node *> queue; 
	std::set<MSG_Node *> records;
	queue.push(&y); records.insert(&y);

	while (!queue.empty()) {
		MSG_Node & next = *(queue.front());
		queue.pop(); VS.erase(&next);

		const MSG_Port & port = y.get_in_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & prev = edge.get_source();
			if (records.count(&prev) == 0
				&& VS.count(&prev) > 0) {
				queue.push(&prev);
				records.insert(&prev);
			}
		}
	}

	return true;
}
bool MSG_Build_Fast::erase_subsumeds(MSG_Node & y, std::set<MSG_Node *> & VS) {
	std::queue<MSG_Node *> queue;
	queue.push(&y); 

	while (!queue.empty()) {
		MSG_Node & next = *(queue.front());
		queue.pop(); 

		const MSG_Port & port = y.get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();
			if (VS.count(&next) > 0) {
				queue.push(&next);
				VS.erase(&next);
			}
		}
	}

	return true;
}
bool MSG_Build_Fast::directSubsumed_randomy(MSG_Node & x, const std::set<MSG_Node *> & base, std::set<MSG_Node *> & DS) {
	std::set<MSG_Node *> VS(base);

	DS.clear();
	while (!VS.empty()) {
		MSG_Node & y = *(*(VS.begin()));
		VS.erase(&y);

		if (subsume(x, y)) {
			DS.insert(&y);
			erase_subsumeds(y, VS);
		}
		else erase_subsuming(y, VS);
	}

	/* purify DS */ return purify_subsumeds(DS);
}

bool MSG_Build_Quick::construct() {
	clustering();
	return linking();
}
bool MSG_Build_Quick::clustering() {
	BitTrieTree trie; ScoreVector * vec;

	clusters.clear();
	while ((vec = producer->produce()) != nullptr) {
		Mutant::ID mid = vec->get_mutant();
		const BitSeq & bits = vec->get_vector();

		BitTrie * leaf = trie.insert_vector(bits);
		if (leaf->get_data() == nullptr) {
			MSG_Node & node = graph.new_node(bits);
			leaf->set_data(&node);
			graph.add_mutant(node, mid);
			clusters.insert(&node);
		}
		else {
			MSG_Node & node = *((MSG_Node *)(leaf->get_data()));
			graph.add_mutant(node, mid);
		}

		consumer->consume(vec);
	}

	return true;
}
bool MSG_Build_Quick::linking() {
	this->build_up(clusters);
	clusters.clear(); return true;
}
bool MSG_Build_Quick::build_up(const std::set<MSG_Node *> & nodes) {
	if (nodes.empty()) return true;		/* do nothing */
	else {
		/* classify nodes for x */
		std::set<MSG_Node *> S, D, N;
		MSG_Node & x = *(*(nodes.begin()));
		this->classify(nodes, x, S, D, N);

		/* recursively build up sub-graph */
		build_up(S); build_up(D); build_up(N);

		/* get the roots and leafs for S, D, N */
		std::set<MSG_Node *> sleafs, sroots;
		std::set<MSG_Node *> dleafs, droots;
		std::set<MSG_Node *> nleafs, nroots;
		derive_roots_leafs(S, sroots, sleafs);
		derive_roots_leafs(D, droots, dleafs);
		derive_roots_leafs(N, nroots, nleafs);

		/* compute direct subsumption */
		std::map<MSG_Node *, std::set<MSG_Node *> *> SN;
		std::map<MSG_Node *, std::set<MSG_Node *> *> ND;
		combine_LR(sleafs, nroots, SN);
		combine_LR(nleafs, droots, ND);

		/* build up direct subsumption */
		// STEP 1. sleafs --> x
		auto sleafs_beg = sleafs.begin();
		auto sleafs_end = sleafs.end();
		while (sleafs_beg != sleafs_end) {
			MSG_Node & src = *(*(sleafs_beg++));
			graph.connect(src, x);
		}
		// STEP 2. x --> droots
		auto droots_beg = droots.begin();
		auto droots_end = droots.end();
		while (droots_beg != droots_end) {
			MSG_Node & trg = *(*(droots_beg++));
			graph.connect(x, trg);
		}
		// STEP 3. S --> N
		auto SN_beg = SN.begin();
		auto SN_end = SN.end();
		while (SN_beg != SN_end) {
			MSG_Node & src = *(SN_beg->first);
			std::set<MSG_Node *> & nexts = *(SN_beg->second);

			auto nexts_beg = nexts.begin();
			auto nexts_end = nexts.end();
			while (nexts_beg != nexts_end) {
				MSG_Node & trg = *(*(nexts_beg++));
				graph.connect(src, trg);
			}

			SN_beg++; delete &nexts;
		}
		// STEP 4. N --> D
		auto ND_beg = ND.begin();
		auto ND_end = ND.end();
		while (ND_beg != ND_end) {
			MSG_Node & src = *(ND_beg->first);
			std::set<MSG_Node *> & nexts = *(ND_beg->second);

			auto nexts_beg = nexts.begin();
			auto nexts_end = nexts.end();
			while (nexts_beg != nexts_end) {
				MSG_Node & trg = *(*(nexts_beg++));
				graph.connect(src, trg);
			}

			ND_beg++; delete &nexts;
		}

		/* end */ return true;
	}
}

bool MSG_Build_Quick::classify(const std::set<MSG_Node *> & C,
	MSG_Node & I, std::set<MSG_Node *> & S,
	std::set<MSG_Node *> & D, std::set<MSG_Node *> & N) {
	S.clear(); D.clear(); N.clear();

	auto beg = C.begin(), end = C.end();
	while (beg != end) {
		MSG_Node & x = *(*(beg++));
		if (&x == &I) continue;
		else if (subsume(I, x)) D.insert(&x);
		else if (subsume(x, I)) S.insert(&x);
		else N.insert(&x);
	}

	return true;
}
bool MSG_Build_Quick::derive_roots_leafs(const std::set<MSG_Node *> & nodes,
	std::set<MSG_Node *> & roots, std::set<MSG_Node *> & leafs) {
	roots.clear(); leafs.clear();

	auto beg = nodes.begin();
	auto end = nodes.end();
	while (beg != end) {
		MSG_Node & x = *(*(beg++));
		if (x.get_in_port().degree() == 0) roots.insert(&x);
		if (x.get_ou_port().degree() == 0) leafs.insert(&x);
	}

	return true;
}
bool MSG_Build_Quick::combine_LR(const std::set<MSG_Node *> & aleafs,
	const std::set<MSG_Node *> & broots,
	std::map<MSG_Node *, std::set<MSG_Node *> *> & ans) {
	/* initialization */
	std::queue<MSG_Node *> squeue;
	std::set<MSG_Node *> svisits;
	auto sbeg = aleafs.begin();
	auto send = aleafs.end();
	while (sbeg != send) {
		squeue.push(*(sbeg));
		svisits.insert(*sbeg);
		sbeg++;
	}

	/* iterate */
	ans.clear(); std::set<MSG_Node *> DS;
	while (!squeue.empty()) {
		/* get next node in A-graph */
		MSG_Node & source = *(squeue.front()); squeue.pop();

		/* compute its direct subsumed nodes */
		derive_direct_subsumed(source, broots, DS);
		update_direct_subsumed(source, DS);

		/* record the answer for source */
		if (!DS.empty())
			ans[&source] = new std::set<MSG_Node *>(DS);

		/* to the upper level */
		const MSG_Port & port = source.get_in_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & prev = edge.get_source();
			if (svisits.count(&prev) == 0) {
				svisits.insert(&prev);
				squeue.push(&prev);
			}
		}
	}

	DS.clear(); return true;
}
bool MSG_Build_Quick::topdown_ready(MSG_Node & next, const std::set<MSG_Node *> & visits) {
	const MSG_Port & port = next.get_in_port();
	for (int i = 0; i < port.degree(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & prev = edge.get_source();
		if (visits.count(&prev) == 0) return false;
	}
	return true;
}
bool MSG_Build_Quick::derive_direct_subsumed(MSG_Node & source,
	const std::set<MSG_Node *> &troots, std::set<MSG_Node *> & DS) {
	/* initialization */
	std::queue<MSG_Node *> tqueue;
	std::set<MSG_Node *> tvisits, records;
	auto tbeg = troots.begin();
	auto tend = troots.end();
	while (tbeg != tend) {
		records.insert(*tbeg);
		tqueue.push(*(tbeg++));
	}

	DS.clear();
	while (!tqueue.empty()) {
		/* get the next node in queue */
		MSG_Node * target = tqueue.front();
		tqueue.pop(); tvisits.insert(target);

		if (subsume(source, *target))
			DS.insert(target);
		/* put the target's "valid" child to queue
		1) not in queue yet;
		2) all the parents have been visited.
		*/
		else {
			const MSG_Port & port = target->get_ou_port();
			for (int i = 0; i < port.degree(); i++) {
				MSG_Edge & edge = port.get_edge(i);
				MSG_Node & next = edge.get_target();
				if (records.count(&next) == 0 &&
					topdown_ready(next, tvisits)) {
					tqueue.push(&next);
					records.insert(&next);
				}
			}
		}
	}

	return true;
}
bool MSG_Build_Quick::update_direct_subsumed(MSG_Node & source, std::set<MSG_Node *> & DS) {
	std::set<MSG_Node *> IS;
	const MSG_Port & port = source.get_ou_port();
	for (int i = 0; i < port.degree() && !DS.empty(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & next = edge.get_target();

		auto beg = DS.begin();
		auto end = DS.end();
		while (beg != end) {
			MSG_Node & check = *(*(beg++));
			if (subsume(next, check)) {
				IS.insert(&check);
			}
		}

		beg = IS.begin(), end = IS.end();
		while (beg != end) DS.erase(*(beg++));
		IS.clear();
	}
	return purify_direct_subsumed(DS);
}
bool MSG_Build_Quick::purify_direct_subsumed(std::set<MSG_Node *> & DS) {
	std::set<MSG_Node *> trash;
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MSG_Node * src = *(beg++);
		if (trash.count(src) > 0) continue;
		else {
			auto beg2 = DS.begin();
			while (beg2 != end) {
				MSG_Node * trg = *(beg2++);
				if (src == trg) continue;
				else if (trash.count(trg) > 0) continue;
				else if (subsume(*src, *trg))
					trash.insert(trg);
			}
		}
	}

	beg = trash.begin(), end = trash.end();
	while (beg != end) DS.erase(*(beg++));
	return true;
}

const std::set<MSG_Pair *> & MSG_Relation::get_related_sources(MSG_Node & node) const {
	if (trg_src.count(node.get_node_id()) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"MSG_Relation::get_rel;ated_sources", 
			"Undefined node: " + std::to_string(node.get_node_id()));
		CErrorConsumer::consume(error); 
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = trg_src.find(node.get_node_id());
		return *(iter->second);
	}
}
const std::set<MSG_Pair *> & MSG_Relation::get_related_targets(MSG_Node & node) const {
	if (src_trg.count(node.get_node_id()) == 0) {
		CError error(CErrorType::InvalidArguments,
			"MSG_Relation::get_rel;ated_sources",
			"Undefined node: " + std::to_string(node.get_node_id()));
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = src_trg.find(node.get_node_id());
		return *(iter->second);
	}
}
void MSG_Relation::build_up() {
	MutantSpace & mspace = source.get_space();
	Mutant::ID n = mspace.number_of_mutants();
	for (Mutant::ID i = 0; i < n; i++) {
		if (source.has_node_of(i) && target.has_node_of(i)) {
			/* get source | target nodes */
			MSG_Node & src = source.get_node_of(i);
			MSG_Node & trg = target.get_node_of(i);
			if (src.get_score_degree() == 0) continue;
			if (trg.get_score_degree() == 0) continue;

			/* compute their key */
			std::string key = "[" + std::to_string(src.get_node_id()) 
				+ "][" + std::to_string(trg.get_node_id()) + "]";
			
			/* construct the maps */
			if (pairs.count(key) == 0) {
				MSG_Pair * pair = new MSG_Pair(src, trg);
				pairs[key] = pair;

				if (src_trg.count(src.get_node_id()) == 0)
					src_trg[src.get_node_id()] = new std::set<MSG_Pair *>();
				if (trg_src.count(trg.get_node_id()) == 0)
					trg_src[trg.get_node_id()] = new std::set<MSG_Pair *>();

				auto iter1 = src_trg.find(src.get_node_id());
				(iter1->second)->insert(pair);
				auto iter2 = trg_src.find(trg.get_node_id());
				(iter2->second)->insert(pair);
			}

			/* insert mutants */
			auto iter = pairs.find(key);
			(iter->second)->add_mutant(i);
		} // end if
	} // end for
}
void MSG_Relation::clear_all() {
	auto beg = pairs.begin();
	auto end = pairs.end();
	while (beg != end)
		delete ((beg++)->second);
	pairs.clear(); 
	src_trg.clear();
	trg_src.clear();
}

void MSG_Tester::gen_tests(const std::set<MSG_Node *> & nodes, TestSet & tests) {
	// initialization
	tests.clear();
	std::set<MSG_Node *> tnodes(nodes);

	while (!tnodes.empty()) {
		// get next node in set 
		MSG_Node & node = *(*(tnodes.begin()));
		tnodes.erase(&node);
		if (node.get_score_degree() == 0) continue;
		const BitSeq & bits = node.get_score_vector();

		// select the next test
		TestCase::ID tid;
		for (tid = 0; tid < bits.bit_number(); tid++) {
			if (bits.get_bit(tid) == BIT_1) {
				tests.add_test(tid); break;
			}
		}

		// remove killed nodes
		if (tid < bits.bit_number())
			eliminate(tid, tnodes);
	}
}
void MSG_Tester::gen_tests(const std::set<Mutant::ID> 
	& mutants, TestSet & tests) {
	// collect requirements
	std::set<MSG_Node *> nodes;
	auto beg = mutants.begin();
	auto end = mutants.end();
	while (beg != end) {
		Mutant::ID mid = *(beg++);
		if (graph->has_node_of(mid)) {
			MSG_Node & node = 
				graph->get_node_of(mid);
			nodes.insert(&node);
		}
	}

	// generate tests
	gen_tests(nodes, tests);
}
void MSG_Tester::eliminate(TestCase::ID tid, std::set<MSG_Node *> & nodes) {
	std::set<MSG_Node *> removed;
	auto beg = nodes.begin();
	auto end = nodes.end();
	while (beg != end) {
		MSG_Node * node = *(beg++);
		const BitSeq & bits = node->get_score_vector();
		if (bits.get_bit(tid) == BIT_1) 
			removed.insert(node);
	}

	beg = removed.begin(), end = removed.end();
	while (beg != end) nodes.erase(*(beg++));
}
bool MSG_Tester::is_killed(const MSG_Node & node, const TestSet & tests) {
	if (node.get_score_degree() == 0) return false;
	else {
		const BitSeq & scores = node.get_score_vector();
		const BitSeq & testes = tests.get_set_vector();
		BitSeq result(scores); result.conjunct(testes);
		return !result.all_zeros();
	}
}
double MSG_Tester::eval_score(const TestSet & tests) {
	size_t total = 0, killed = 0;
	long n = graph->size(), k;
	for (k = 0; k < n; k++) {
		MSG_Node & node = graph->get_node(k);
		if (node.get_score_degree() == 0) continue;

		size_t num = node.get_mutants().number_of_mutants();
		total += num; 
		if (is_killed(node, tests)) killed += num;
	}
	return ((double)killed) / ((double)total);
}
void MSG_Tester::collect_roots(std::set<MSG_Node *> &roots) {
	roots.clear();

	long k, n = graph->size();
	MSG_Node * eqs = nullptr;
	for (k = 0; k < n; k++) {
		MSG_Node & node = graph->get_node(k);
		if (node.get_score_degree() == 0) {
			eqs = &node; break;
		}
		else if (node.get_in_port().degree() == 0) {
			roots.insert(&node);
		}
	}

	if (eqs != nullptr) {
		roots.clear();
		const MSG_Port & port = eqs->get_ou_port();
		for (int k = 0; k < port.degree(); k++) {
			MSG_Edge & edge = port.get_edge(k);
			MSG_Node & next = edge.get_target();
			roots.insert(&next);
		}
	}
}
double MSG_Tester::eval_dom_score(const TestSet & tests) {
	std::set<MSG_Node *> roots;
	collect_roots(roots);

	size_t killed = 0;
	auto beg = roots.begin();
	auto end = roots.end();
	while (beg != end) {
		MSG_Node & node = *(*(beg++));
		if (is_killed(node, tests)) killed++;
	}

	return ((double)killed) / ((double)roots.size());
}




