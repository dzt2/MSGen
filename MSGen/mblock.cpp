#include "mblock.h"

MutBlock::MutBlock(MutBlockGraph & g, ID id, const BitSeq & covvec)
	: bid(id), coverage(covvec), graph(g), mutants(), local_graph(g.get_mutant_space()), ports() {}
BlockPort & MutBlock::get_port_of(MutBlock & y) const {
	if (ports.count(&y) == 0) {
		CError error(CErrorType::InvalidArguments, "MutBlock::get_port_of", "Undefined port: "
			+ std::to_string(bid) + " --> " + std::to_string(y.get_block_id()));
		CErrorConsumer::consume(error);
	}

	auto iter = ports.find(&y);
	return *(iter->second);
}
void MutBlock::link_to_block(MutBlock & y) {
	if (ports.count(&y) == 0 && (&y) != this) 
		ports[&y] = new BlockPort(*this, y);
}
void MutBlock::dislink_from_block(MutBlock & y) {
	if (ports.count(&y) == 0) {
		CError error(CErrorType::InvalidArguments, "MutBlock::dislink_from_block", "Undefined port: "
			+ std::to_string(bid) + " --> " + std::to_string(y.get_block_id()));
		CErrorConsumer::consume(error);
	}
	else {
		auto iter = ports.find(&y);
		delete iter->second;
		ports.erase(&y);
	}
}
void MutBlock::dislink_from_all() {
	auto beg = ports.begin();
	auto end = ports.end();
	while (beg != end) 
		delete (beg++)->second;
	ports.clear();
}

const std::list<MuSubsume> & BlockPort::get_direct_subsumption(MSGVertex & x) const {
	if (edges.count(&x) == 0) {
		CError error(CErrorType::InvalidArguments, "BlockPort::get_direct_subsumption", "Unlinked node: " + x.get_id());
		CErrorConsumer::consume(error);
	}
	
	auto iter = edges.find(&x); return *(iter->second);
}
void BlockPort::clear() {
	auto beg = edges.begin();
	auto end = edges.end();
	while (beg != end) 
		delete (beg++)->second;
	
	edges.clear();
	valid_nodes.clear();

}
void BlockPort::init() {
	clear();

	const BitSeq & trg_coverage = target.get_coverage();
	MSGraph & msg = source.get_local_graph();
	MSGVertex::ID vid = 0, vnum = msg.number_of_vertices();
	while (vid < vnum) {
		MSGVertex & x = msg.get_vertex(vid++);
		const BitSeq & score_vector = x.get_feature()->get_vector();
		if (score_vector.subsume(trg_coverage)) valid_nodes.insert(&x);
	}
}
void BlockPort::link(MSGVertex & x, MSGVertex & y) {
	if (valid_nodes.count(&x) == 0) {
		CError error(CErrorType::InvalidArguments, "BlockPort::link", 
			"Invalid source node: " + std::to_string(x.get_id()));
		CErrorConsumer::consume(error);
	}
	else {
		std::list<MuSubsume> * DS;
		if (edges.count(&x) == 0) {
			DS = new std::list<MuSubsume>();
			edges[&x] = DS;
		}
		else {
			auto iter = edges.find(&x);
			DS = iter->second;
		}

		MuSubsume edge(x, y);
		DS->push_front(edge);
	}
}
void BlockPort::update() {
	valid_nodes.clear();

	auto beg = edges.begin();
	auto end = edges.end();
	while (beg != end)
		valid_nodes.insert((beg++)->first);
}

MutBlock & MutBlockGraph::get_block(MutBlock::ID bid) const {
	if (bid >= blocks.size()) {
		CError error(CErrorType::InvalidArguments,
			"MutBlockGraph::get_block",
			"Undefined block: " + std::to_string(bid));
		CErrorConsumer::consume(error);
	}

	return *(blocks[bid]);
}
MutBlock & MutBlockGraph::get_block_of_mutant(Mutant::ID mid) const {
	if (index.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"MutBlockGraph::get_block_of_mutants", 
			"Undefined mutant: " + std::to_string(mid));
		CErrorConsumer::consume(error);
	}

	auto iter = index.find(mid);
	return *(iter->second);
}
MutBlock * MutBlockGraph::new_block(const BitSeq & covvec) {
	MutBlock * block = new MutBlock(*this, blocks.size(), covvec);
	blocks.push_back(block); return block;
}
void MutBlockGraph::add_mutant(MutBlock * block, Mutant::ID mid) {
	if (index.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments,
			"MutBlockGraph::add_mutant",
			"Undefined mutant: " + std::to_string(mid));
		CErrorConsumer::consume(error);
	}
	else { index[mid] = block; block->add_mutant(mid); }
}
void MutBlockGraph::clear() {
	auto beg = blocks.begin();
	auto end = blocks.end();
	while (beg != end) delete *(beg++);
	index.clear(); blocks.clear();
}

void MutBlockBuilder::open(MutBlockGraph & bgraph) {
	close(); 
	graph = &bgraph;
	trie = new BitTrieTree();
}
void MutBlockBuilder::close() {
	if (graph != nullptr) {
		delete trie;
		trie = nullptr;
		graph = nullptr;
	}
}
void MutBlockBuilder::build(const CoverageVector & covvec) {
	/* getters */
	Mutant::ID mid = covvec.get_mutant();
	const BitSeq & coverage = covvec.get_coverage();

	if (coverage.all_zeros()) return;	/* filter */

	/* get the block for the coverage */
	BitTrie * leaf = trie->insert_vector(coverage);
	MutBlock * block;
	if (leaf->get_data() == nullptr) {
		block = graph->new_block(coverage);
		leaf->set_data(block);
	}
	else block = (MutBlock *)(leaf->get_data());

	/* add mutant into the block */
	graph->add_mutant(block, mid);
}
void MutBlockBuilder::build_blocks(MutBlockGraph & bgraph, 
	CoverageProducer & producer, CoverageConsumer & consumer) {
	CoverageVector * covvec;

	this->open(bgraph);
	while ((covvec = producer.produce()) != nullptr) 
		this->build(*covvec);
	this->close();
}

void LocalMSGBuilder::open(MutBlockGraph & bgraph) {
	close();
	graph = &bgraph;
}
void LocalMSGBuilder::close() { 
	auto beg = builders.begin();
	auto end = builders.end();
	while (beg != end) delete (beg++)->second;
	builders.clear(); graph = nullptr;
}
void LocalMSGBuilder::add(const ScoreVector & score_vector) {
	/* get mutant */
	Mutant::ID mid = score_vector.get_mutant();
	if (score_vector.get_vector().all_zeros()) return;

	/* get builder */
	MSGBuilder * builder;
	MutBlock & block = graph->get_block_of_mutant(mid);
	if (builders.count(&block) == 0) {
		builder = new MSGBuilder();
		builder->open(block.get_local_graph());
		builders[&block] = builder;
	}
	else {
		auto iter = builders.find(&block);
		builder = iter->second;
	}

	/* cluster */
	builder->add_score_vector(score_vector);
}
void LocalMSGBuilder::end() {
	auto beg = builders.begin();
	auto end = builders.end();
	while (beg != end) {
		MSGBuilder & builder = *((beg++)->second);
		builder.end_score_vectors();
	}
}
void LocalMSGBuilder::construct_local_graph(MutBlockGraph & bgraph, 
	ScoreProducer & producer, ScoreConsumer & consumer) {
	ScoreVector * score_vector;

	this->open(bgraph);
	while ((score_vector = producer.produce()) != nullptr) {
		this->add(*score_vector);
	}
	this->end(); this->close();
}

void MutBlockConnect::connect(MutBlockGraph & bgraph) {
	MutBlock::ID bid = 0, bnum = bgraph.number_of_blocks();

	this->open(bgraph);
	this->init_connection();
	while (bid < bnum) {
		MutBlock & block = bgraph.get_block(bid++);
		this->build_up_ports(block);
	}
	this->close();
}
void MutBlockConnect::open(MutBlockGraph & bgraph) {
	close(); graph = &bgraph;
}
void MutBlockConnect::init_connection() {
	MutBlock::ID i, j, num = graph->number_of_blocks();
	for (i = 0; i < num; i++) {
		MutBlock & bi = graph->get_block(i);
		const BitSeq & covi = bi.get_coverage();
		bi.dislink_from_all();
		for (j = 0; j < num; j++) {
			if (i == j) continue;
			MutBlock & bj = graph->get_block(j);
			BitSeq  covj = bj.get_coverage();
			covj.conjunct(covi);
			if (!covj.all_zeros()) 
				bi.link_to_block(bj);
		}
	}
}
void MutBlockConnect::build_up_ports(MutBlock & block) {
	const std::map<MutBlock *, BlockPort *> & ports = block.get_ports();
	std::set<MutBlock *> invalid_targets;

	auto beg = ports.begin(), end = ports.end();
	while (beg != end) {
		/* get the target and its port to it */
		MutBlock * trg = beg->first;
		BlockPort & port = *((beg++)->second);
		/* compute the subsumption for port */
		this->build_up_port(port);
		/* if port is empty, record target as invalid */
		if (port.get_valid_nodes().empty()) 
			invalid_targets.insert(trg);
	}

	/* dislink the block from invalid targets and remove their ports */
	auto tbeg = invalid_targets.begin();
	auto tend = invalid_targets.end();
	while (tbeg != tend) {
		MutBlock & trg = *(*(tbeg++));
		block.dislink_from_block(trg);
	}
}
void MutBlockConnect::build_up_port(BlockPort & port) {
	MutPortComputer linker(port);
	linker.compute();
}
void MutBlockConnect::close() { graph = nullptr; }

MutPortComputer::~MutPortComputer() {
	while (!questions.empty()) questions.pop();

	auto beg = solutions.begin();
	auto end = solutions.end();
	while (beg != end) delete (beg++)->second;
	solutions.clear(); leafs.clear();
}
void MutPortComputer::compute() {
	this->open_questions();
	while (!questions.empty()) {
		/* get next node */
		MSGVertex * x = questions.front();
		questions.pop();

		/* invalid access */
		if (solutions.count(x) > 0) {
			CError error(CErrorType::Runtime, "MutPortComputer::compute", 
				"Duplicated question for msg-vertex (" + std::to_string(x->get_id()) + ")");
			CErrorConsumer::consume(error);
		}
		/* solving the question */
		else {
			std::set<MSGVertex *> * DS = new std::set<MSGVertex *>();
			this->initial_subsumption(x, *DS);
			this->compute_subsumption(x, *DS);
			this->connect_subsumption(x, *DS);
			solutions[x] = DS;
		}
	}
	this->close_questions();
}
void MutPortComputer::open_questions() {
	if (!questions.empty() || !solutions.empty()) {
		CError error(CErrorType::Runtime, 
			"MutPortComputer::open_questions", 
			"Invalid access to open questions");
		CErrorConsumer::consume(error);
	}
	else {
		leafs.clear();	/* initial leafs for valid nodes */
		std::queue<MSGVertex *> qlist;	/* to similate the question list */
		std::set<MSGVertex *> visits;	/* those have been "solved" */

		/* insert leafs in valid node-graph to the queue */
		const std::set<MSGVertex *> & valid_nodes = port.get_valid_nodes();
		auto beg = valid_nodes.begin(), end = valid_nodes.end();
		while (beg != end) {
			MSGVertex * node = *(beg++);
			if (is_leaf(node, valid_nodes)) {
				qlist.push(node); 
				visits.insert(node);
				leafs.insert(node);
			}
		}

		/* update the question queue */
		while (!qlist.empty()) {
			MSGVertex * x = qlist.front(); 
			qlist.pop(); questions.push(x);

			const std::list<MuSubsume> & edges 
				= x->get_in_port().get_edges();
			auto ebeg = edges.begin();
			auto eend = edges.end();
			while (ebeg != eend) {
				const MuSubsume & edge = *(ebeg++);
				MSGVertex * parent = (MSGVertex *)(&(edge.get_source()));
				if (visits.count(parent) > 0) continue;
				else if (is_solvable(parent, visits, valid_nodes)) {
					qlist.push(parent); visits.insert(parent);
				}
			}
		} /* end while: qlist.empty() */

		/* clear solutions */ solutions.clear();
	}
}
void MutPortComputer::initial_subsumption(MSGVertex * x, std::set<MSGVertex *> & DS) {
	/* initial seed for leaf is the subsumed leaf of another MSG in target block */
	/*if (leafs.count(x) > 0) {
		const std::set<MSGVertex *> & msg_leafs
			= port.get_target_block().get_local_graph().get_leafs();
		const BitSeq & xvec = x->get_feature()->get_vector();

		auto beg = msg_leafs.begin(), end = msg_leafs.end();
		while (beg != end) {
			MSGVertex * y = *(beg++);
			const BitSeq & yvec = y->get_feature()->get_vector();
			if (xvec.subsume(yvec)) DS.insert(y);
		}
	}*/
	/* initial seed for non-leaf is the solution by its valid targets */
	/*else {
		const std::list<MuSubsume> & edges = x->get_ou_port().get_edges();
		const std::set<MSGVertex *> & valid_nodes = port.get_valid_nodes();

		auto edge_beg = edges.begin(), edge_end = edges.end();
		while (edge_beg != edge_end) {
			const MuSubsume & edge = *(edge_beg++);
			MSGVertex * y = (MSGVertex *)(&(edge.get_target()));
			if (valid_nodes.count(y) == 0) continue;

			auto iter = solutions.find(y);
			std::set<MSGVertex *> * y_DS = iter->second;
			this->add_all(*y_DS, DS);
		}
	}*/

	const std::set<MSGVertex *> & msg_leafs
		= port.get_target_block().get_local_graph().get_leafs();
	const BitSeq & xvec = x->get_feature()->get_vector();

	auto beg = msg_leafs.begin(), end = msg_leafs.end();
	while (beg != end) {
		MSGVertex * y = *(beg++);
		const BitSeq & yvec = y->get_feature()->get_vector();
		if (xvec.subsume(yvec)) DS.insert(y);
	}
}
void MutPortComputer::compute_subsumption(MSGVertex * x, std::set<MSGVertex *> & DS) {
	/* declarations */
	std::queue<MSGVertex *> subsuming_queue;
	std::set<MSGVertex *> subsumes, nonsubsumes;

	/* initialization */
	auto DS_beg = DS.begin(), DS_end = DS.end();
	while (DS_beg != DS_end) {
		MSGVertex * node = *(DS_beg++);
		subsuming_queue.push(node);
		subsumes.insert(node);
	}
	DS.clear(); bool direct_subsume;

	/* compute by the subsuming queue */
	while (!subsuming_queue.empty()) {
		/* get next subsumed node by x */
		MSGVertex * y = subsuming_queue.front();
		subsuming_queue.pop(); direct_subsume = true;

		/* get nodes subsuming y in target block */
		const std::list<MuSubsume> & edges 
			= y->get_in_port().get_edges();
		auto edge_beg = edges.begin();
		auto edge_end = edges.end();
		while (edge_beg != edge_end) {
			/* get next parent of the y in target block */
			const MuSubsume & edge = *(edge_beg++); bool is_subsume;
			MSGVertex * parent = (MSGVertex *)(&(edge.get_source()));

			/* compute subsumption from x to parent */
			if (subsumes.count(parent) > 0)
				is_subsume = true;
			else if (nonsubsumes.count(parent) > 0)
				is_subsume = false;
			else is_subsume = subsume(*x, *parent);

			/* if x subsumes y's parent:
				1) x does not directly subsume y;
				2) push y's parent to the queue for further analysis
			*/
			if (is_subsume) {
				direct_subsume = false;
				if (subsumes.count(parent) == 0) {
					subsumes.insert(parent);
					subsuming_queue.push(parent);
				}
			}
			/* otherwise:
				1) x does not subsume y's parent;
				2) its parent is recorded;
				3) its parent will not be inserted to the queue
			*/
			else nonsubsumes.insert(parent);
		} /* end while: edge_beg --> edge_end */

		/* directly subsumed */
		if (direct_subsume) DS.insert(y);
	} /* end while: */

	/* remove all those subsumed by its children */
	const std::list<MuSubsume> & out_edges = x->get_ou_port().get_edges();
	auto out_beg = out_edges.begin(), out_end = out_edges.end();
	while (out_beg != out_end) {
		const MuSubsume & edge = *(out_beg++);
		MSGVertex * child = (MSGVertex *)(&(edge.get_target()));
		if (solutions.count(child) > 0) {
			auto iter = solutions.find(child);
			std::set<MSGVertex *> * SDS = iter->second;
			this->sub_set(DS, *SDS);
		}
	}

	/* return */ return;
}
void MutPortComputer::connect_subsumption(MSGVertex * x, std::set<MSGVertex *> & DS) {
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MSGVertex * y = *(beg++);
		port.link(*x, *y);
	}
}
void MutPortComputer::close_questions() { port.update(); }
bool MutPortComputer::is_leaf(MSGVertex * x, const std::set<MSGVertex *> & valid_nodes) {
	if (valid_nodes.count(x) > 0) {
		const std::list<MuSubsume> & edges = x->get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			const MuSubsume & edge = *(beg++);
			MSGVertex * y = (MSGVertex *)(&(edge.get_target()));
			if (valid_nodes.count(y) > 0) return false;
		}
		return true;
	}
	else return false;
}
bool MutPortComputer::is_solvable(MSGVertex * x, const std::set<MSGVertex *> & solutions, const std::set<MSGVertex *> & valid_nodes) {
	if (valid_nodes.count(x) > 0) {
		const std::list<MuSubsume> & edges = x->get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			const MuSubsume & edge = *(beg++);
			MSGVertex * y = (MSGVertex *)(&(edge.get_target()));

			/* invalid target is not considered */
			if (valid_nodes.count(y) == 0) continue;
			/* valid target is not solved, x is unsolvable */
			else if (solutions.count(y) == 0) return false;
		}
		return true;
	}
	else return false;
}
void MutPortComputer::add_all(std::set<MSGVertex *> & src, std::set<MSGVertex *> & trg) {
	auto beg = src.begin();
	auto end = src.end();
	while (beg != end) 
		trg.insert(*(beg++));
}
void MutPortComputer::sub_set(std::set<MSGVertex *> & a, std::set<MSGVertex *> & b) {
	auto beg = b.begin(), end = b.end();
	while (beg != end) {
		MSGVertex * vex = *(beg++);
		if (a.count(vex) > 0) a.erase(vex);
	}
}
bool MutPortComputer::subsume(MSGVertex & x, MSGVertex & y) {
	const BitSeq & xv = x.get_feature()->get_vector();
	const BitSeq & yv = y.get_feature()->get_vector();
	return xv.subsume(yv);
}

/* calculate the number of edges in local graph */
unsigned int number_of_edges(const MSGraph & graph) {
	unsigned int edges = 0;
	MSGVertex::ID vid, vnum = graph.number_of_vertices();
	for (vid = 0; vid < vnum; vid++) {
		MSGVertex & vex = graph.get_vertex(vid);
		edges += vex.get_ou_degree();
	}
	return edges;
}
/* print the blocks to output stream */
static void print_blocks(MutBlockGraph & bgraph, std::ostream & out) {
	/* iterate each block to print their information */
	MutBlock::ID bid = 0, bnum = bgraph.number_of_blocks();

	/* print blocks */
	while (bid < bnum) {
		MutBlock & block = bgraph.get_block(bid++);
		out << "  Block[" << block.get_block_id() << "] = { M: ";
		out << block.get_mutants().size() << "; C: "
			<< block.get_local_graph().number_of_vertices()
			<< "; \tB: " << block.get_coverage().to_string() << " }\n";

		const MSGraph & graph = block.get_local_graph();
		MSGVertex::ID vid = 0, vnum = graph.number_of_vertices();
		out << "  /================== Vertex =================/\n";
		while (vid < vnum) {
			MSGVertex & vex = graph.get_vertex(vid);
			out << "\tV[" << vid << "] = { M: " << vex.get_cluster().number_of_mutants();
			out << "; \tB: " << vex.get_feature()->get_vector().to_string() << " }\n";
			vid++;
		}

		out << "\n  /================== Edges =================/\n";
		for (vid = 0; vid < vnum; vid++) {
			MSGVertex & vex = graph.get_vertex(vid);
			const std::list<MuSubsume> & edges = vex.get_ou_port().get_edges();
			auto ebeg = edges.begin(), eend = edges.end();
			while (ebeg != eend) {
				const MuSubsume & edge = *(ebeg++);
				out << "\tSubsume { " << edge.get_source().get_id();
				out << ", " << edge.get_target().get_id() << " }\n";
			}
		}
		out << "\n\n";
	}

	/* print ports */
	for (bid = 0; bid < bnum; bid++) {
		/* get the ports of the ith block */
		const std::map<MutBlock *, BlockPort *> & ports
			= bgraph.get_block(bid).get_ports();
		auto beg = ports.begin(), end = ports.end();

		/* iterate ports in one block */
		while (beg != end) {
			/* get next port */
			BlockPort & port = *((beg++)->second);
			out << "  Port [ " << port.get_source_block().get_block_id();
			out << ", " << port.get_target_block().get_block_id() << " ]\n";

			const std::set<MSGVertex *> & vertices = port.get_valid_nodes();
			auto vbeg = vertices.begin(), vend = vertices.end();
			while (vbeg != vend) {
				/* get next source node */
				MSGVertex & vertex = *(*(vbeg++));
				const std::list<MuSubsume> & inter_edges 
					= port.get_direct_subsumption(vertex);
				out << "Vex[" << vertex.get_id() << "]: ";

				/* get the edges from source to nodes in target block */
				auto inter_beg = inter_edges.begin();
				auto inter_end = inter_edges.end();
				while (inter_beg != inter_end) {
					const MuSubsume & inter_edge = *(inter_beg++);
					out << inter_edge.get_target().get_id() << "; ";
				}
				out << "\n";
			}

			out << "\n";
		}
	}

	out << std::endl;
	/* return */ return;
}
/* statistic analysis on output stream */
static void analysis_blocks(MutBlockGraph & bgraph, std::ostream & out) {
	out << "Total Statistics on Mutation Blocks\nBlock\t#Mutants\t#Clusters\t#Subsume\n";

	MutBlock::ID bid = 0, bnum = bgraph.number_of_blocks();
	while (bid < bnum) {
		MutBlock & block = bgraph.get_block(bid++);
		const MSGraph & graph = block.get_local_graph();

		out << block.get_block_id() << '\t';
		out << block.get_mutants().size() << "\t";
		out << graph.number_of_vertices() << "\t";
		out << number_of_edges(graph) << "\n";
	}
	out << std::endl;
}
/* statistical analysis on bridges between blocks */
static void analysis_ports(MutBlockGraph & bgraph, std::ostream & out) {
	MutBlock::ID bid = 0, bnum = bgraph.number_of_blocks();

	out << "Total Statistics on Mutation Bridges\nSource\tTarget\tSubsume?\t#Src_Vertex\t#Sub_Vertex\t#Edges\t#Equivalent\t#Strict\n";
	while (bid < bnum) {
		/* get the ports of the next block in the graph */
		MutBlock & src = bgraph.get_block(bid++);
		const std::map<MutBlock *, BlockPort *> & ports = src.get_ports();

		/* iterate all ports in src */
		auto beg = ports.begin(), end = ports.end();
		while (beg != end) {
			BlockPort & port = *((beg++)->second);	/* get next port in the src */

			/* summary */
			out << port.get_source_block().get_block_id() << "\t";
			out << port.get_target_block().get_block_id() << "\t";
			out << port.get_source_block().get_coverage().subsume(port.get_target_block().get_coverage()) << "\t";
			out << port.get_valid_nodes().size() << "\t";

			/* analysis on edges */
			const std::set<MSGVertex *> & vertices = port.get_valid_nodes();
			auto vbeg = vertices.begin(), vend = vertices.end();
			size_t eq = 0, st = 0, sv = 0;
			while (vbeg != vend) {
				MSGVertex & vex = *(*(vbeg++));
				const std::list<MuSubsume> & edges = port.get_direct_subsumption(vex);
				auto ebeg = edges.begin(), eend = edges.end();
				while (ebeg != eend) {
					const MuSubsume & edge = *(ebeg++);
					const BitSeq & svec = edge.get_source().get_feature()->get_vector();
					const BitSeq & tvec = edge.get_target().get_feature()->get_vector();
					if (svec.equals(tvec)) eq++;
					else st++;
				}
				sv++;
			}	/* end while: vbeg != vend */
			out << sv << "\t" << (st + eq) << "\t" << eq << "\t" << st << "\n";
		} /* end while: beg != end */

	} /* end while: blocks iteration */
	out << std::endl;
}

/* validate */
static void validate_equivalence(MutBlockGraph & bgraph) {
	MutBlock::ID bid = 0, bnum = bgraph.number_of_blocks();
	while (bid < bnum) {
		/* get next block for validation */
		MutBlock & block = bgraph.get_block(bid++);
		const std::map<MutBlock *, BlockPort *> & ports = block.get_ports();
		std::cerr << "Block #" << bid - 1 << "\n";

		auto beg = ports.begin(), end = ports.end();
		while (beg != end) {
			/* get next port */
			BlockPort & port = *((beg++)->second); 
			MutBlock::ID tid = port.get_target_block().get_block_id();
			const std::set<MSGVertex *> & valid_nodes = port.get_valid_nodes();

			/* validate its equivalence */
			auto node_beg = valid_nodes.begin();
			auto node_end = valid_nodes.end();
			while (node_beg != node_end) {
				const std::list<MuSubsume> & edges = 
					port.get_direct_subsumption(*(*(node_beg++)));

				auto edge_beg = edges.begin();
				auto edge_end = edges.end();
				while (edge_beg != edge_end) {
					const MuSubsume & edge = *(edge_beg++);
					const MSGVertex & x = edge.get_source();
					const MSGVertex & y = edge.get_target();
					const BitSeq & xvec = x.get_feature()->get_vector();
					const BitSeq & yvec = y.get_feature()->get_vector();

					if (xvec.equals(yvec)) {
						std::cerr << "\tB" << bid - 1 << "[" << x.get_id() << "] --> B" << tid << "[" << y.get_id() << "]\n";
					}
				}
			}
		}
		std::cerr << std::endl;
	}
	
}

int main() { 
	// initialization
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string prname = "prime";
	File & root = *(new File(prefix + prname)); TestType ttype = TestType::general;

	// create code-project, mutant-project, test-project
	CProgram & program = *(new CProgram(root));
	CTest & ctest = *(new CTest(ttype, root, program.get_exec()));
	CMutant & cmutant = *(new CMutant(root, program.get_source()));

	// load tests
	ctest.load(); const TestSpace & tspace = ctest.get_space();
	std::cout << "Loading test cases: " << tspace.number_of_tests() << std::endl;

	// load mutations 
	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		const CodeFile & cfile = *(*(cfile_beg++));
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		cmutant.load_mutants_for(mspace, true);
		std::cout << "Load " << mspace.number_of_mutants() <<
			" mutants for: " << cfile.get_file().get_path() << "\n" << std::endl;
	}

	// get the trace 
	CTrace & ctrace = *(new CTrace(root, cspace, tspace));
	CoverageSpace & covspace = ctrace.get_space();

	// score
	CScore & cscore = *(new CScore(root, cmutant, ctest));

	// builders 
	MutBlockBuilder block_builder;
	LocalMSGBuilder graph_builder;
	MutBlockConnect block_connect;

	// get coverage vector
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
	cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		// get next file and load its text 
		CodeFile & cfile = *(*(cfile_beg++)); cspace.load(cfile);

		// get mutations for code-file
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		MutantSet & mutants = *(mspace.create_set());
		mutants.complement();

		// load coverage 
		covspace.add_file_coverage(cfile, tests);
		ctrace.load_coverage(cfile);
		std::cout << "Loading coverage for \"" << cfile.get_file().get_path() << "\"\n";

		// get coverage producer and consumer
		FileCoverage & fcov = covspace.get_file_coverage(cfile);
		CoverageProducer cproducer(mspace, fcov);
		CoverageConsumer ccomsumer;

		// get score producer and consumer
		ScoreSource & score_src = cscore.get_source(cfile);
		ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
		ScoreProducer sproducer(score_func); ScoreConsumer sconsumer(score_func);

		// create mutblock and builder
		MutBlockGraph & blocks = *(new MutBlockGraph(mspace));

		// building blocks
		block_builder.build_blocks(blocks, cproducer, ccomsumer);
		graph_builder.construct_local_graph(blocks, sproducer, sconsumer);
		block_connect.connect(blocks);

		// print block information 
		std::ofstream fout(prefix + prname + "/bgraph.txt");
		print_blocks(blocks, fout); fout.close();
		std::ofstream bout(prefix + prname + "/blocks.txt");
		analysis_blocks(blocks, bout); bout.close();
		std::ofstream lout(prefix + prname + "/bridges.txt");
		analysis_ports(blocks, lout); lout.close();

		// validation 
		validate_equivalence(blocks);

		// delete resources
		delete &blocks;
	}

	// delete resources
	delete & cscore; delete &ctrace;
	delete &cmutant; delete & ctest;
	delete &program; delete &  root;
	// end to return 
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}