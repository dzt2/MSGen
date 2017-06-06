#include "cblock.h"

size_t MutBlockBridge::get_degree_of(MSGVertex & vex) const {
	/* validate the vertex for input */
	if (vertices.count(&vex) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"MutBlockBridge::get_degree_of", 
			"Invalid vex (" + std::to_string(vex.get_id()) + ")");
		CErrorConsumer::consume(error);
	}

	/* compute its degree */
	if (edges.count(&vex) == 0) return 0;
	else {
		auto iter = edges.find(&vex);
		return (iter->second)->size();
	}
}
const std::list<MuSubsume> & MutBlockBridge::get_edges_of(MSGVertex & vex) const {
	if (edges.count(&vex) == 0) {
		CError error(CErrorType::InvalidArguments,
			"MutBlockBridge::get_edges_of",
			"Invalid vex (" + std::to_string(vex.get_id()) + ")");
		CErrorConsumer::consume(error);
	}

	auto iter = edges.find(&vex);
	return *(iter->second);
}
void MutBlockBridge::init() {
	clear();	/* clear the original edges and nodes */

	const BitSeq & trg_coverage = target.get_coverage();
	const MSGraph & src_graph = source.get_local_graph();
	Mutant::ID vid, vnum = src_graph.number_of_vertices();
	for (vid = 0; vid < vnum; vid++) {
		MSGVertex & vex = src_graph.get_vertex(vid);
		if (vex.get_feature() == nullptr) continue;

		const BitSeq & score_vector = 
			vex.get_feature()->get_vector();
		if (score_vector.subsume(trg_coverage))
			vertices.insert(&vex);
	}
}
void MutBlockBridge::clear() {
	/* delete all edges */
	auto beg = edges.begin();
	auto end = edges.end();
	while (beg != end) {
		auto list = (beg++)->second;
		list->clear(); delete list;
	}

	/* clear candidates */
	vertices.clear();
	edges.clear();
}
void MutBlockBridge::links(MSGVertex & src, MSGVertex & trg) {
	if (vertices.count(&src) == 0) {
		CError error(CErrorType::InvalidArguments, "MutBlockBridge:links", 
			"Invalid src (" + std::to_string(src.get_id()) + ")");
		CErrorConsumer::consume(error);
	}
	else {
		/* get the edges from the specified node */
		std::list<MuSubsume> * targets;
		if (edges.count(&src) == 0) {
			targets = new std::list<MuSubsume>();
			edges[&src] = targets;
		}
		else {
			auto iter = edges.find(&src);
			targets = iter->second;
		}

		/* create edge and put into list */
		MuSubsume edge(src, trg);
		targets->push_front(edge);
	}
}

MutBlock & MutBlockGraph::get_block_of(Mutant::ID mid) const {
	if (index.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, "MutBlockGraph::get_block_of",
			"Undefined mutant (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error);
	}
	
	auto iter = index.find(mid);
	return *(iter->second);
}
MutBlock * MutBlockGraph::new_block(const BitSeq & coverage) {
	MutBlock * block = new MutBlock(coverage, mut_space);
	blocks.insert(block); return block;
}
void MutBlockGraph::add_index(Mutant::ID mid, MutBlock * block) {
	if (index.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments, "MutBlockGraph::add_index",
			"Duplicated mutant for (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error);
	}
	else index[mid] = block;
}
MutBlockBridge * MutBlockGraph::new_bridge(MutBlock & src, MutBlock & trg) {
	/* validation */
	if (blocks.count(&src) == 0) {
		CError error(CErrorType::InvalidArguments, "MutBlockGraph::new_bridge",
			"Undefined block for (src)");
		CErrorConsumer::consume(error);
	}
	else if (blocks.count(&trg) == 0) {
		CError error(CErrorType::InvalidArguments, "MutBlockGraph::new_bridge",
			"Undefined block for (trg)");
		CErrorConsumer::consume(error);
	}
	
	/* validate intersection */
	BitSeq inter = src.get_coverage();
	inter.conjunct(trg.get_coverage());
	if (inter.all_zeros()) {
		CError error(CErrorType::Runtime, "MutBlockGraph::new_bridge",
			"Invalid block pair: (src * trg = \emptyset)");
		CErrorConsumer::consume(error);
	}
	
	/* create bridge and put into list */
	MutBlockBridge * bridge = new MutBlockBridge(src, trg);
	bridges.insert(bridge); return bridge;
}
void MutBlockGraph::clear_bridges() {
	auto brg_beg = bridges.begin();
	auto brg_end = bridges.end();
	while (brg_beg != brg_end)
		delete *(brg_beg++);
	bridges.clear();
}
void MutBlockGraph::clear_all() {
	/* clear bridges between blocks */
	this->clear_bridges();

	/* index are cleared */ index.clear();

	/* blocks are cleared */
	auto blk_beg = blocks.begin();
	auto blk_end = blocks.end();
	while (blk_beg != blk_end)
		delete *(blk_beg++);
	blocks.clear();
}

void MutBlockBuilder::open(MutBlockGraph & graph) {
	close();
	
	bgraph = &graph; 
	graph.clear_all();
	trie = new BitTrieTree();
}
void MutBlockBuilder::close() {
	if (trie != nullptr) {
		delete trie; trie = nullptr;
	}
	bgraph = nullptr;
}
void MutBlockBuilder::insert(const CoverageVector & cvector) {
	if (bgraph != nullptr) {
		/* get the leaf where block is inserted */
		BitTrie * leaf = trie->insert_vector(cvector.get_coverage());
		
		/* get the target block */
		MutBlock * block;
		if (leaf->get_data() == nullptr) {
			block = bgraph->new_block(
				cvector.get_coverage());
			leaf->set_data(block);
		}
		else block = (MutBlock *)(leaf->get_data());

		/* insert mutation */
		block->add_mutant(cvector.get_mutant());
		bgraph->add_index(cvector.get_mutant(), block);
	}
	else {
		CError error(CErrorType::Runtime, 
			"MutBlockBuilder::insert", "Invalid access: not opened");
		CErrorConsumer::consume(error);
	}
}
void MutBlockBuilder::build_blocks(MutBlockGraph & graph, 
	CoverageProducer & producer, CoverageConsumer & consumer) {
	CoverageVector * vector;

	this->open(graph);
	while ((vector = producer.produce()) != nullptr) {
		this->insert(*vector);
		consumer.consume(vector);
	}
	this->close();
}

void LocalMSGBuilder::open(MutBlockGraph & bgraph) {
	/* initialization */
	close(); graph = &bgraph;
	
	/* clear graphs */
	const std::set<MutBlock *> & blocks = bgraph.get_blocks();
	auto beg = blocks.begin(), end = blocks.end();
	while (beg != end) {
		MutBlock & block = *(*(beg++));
		MSGraph & local_graph = block.get_graph_for_modify();
		local_graph.clear();
	}
}
void LocalMSGBuilder::close() {
	graph = nullptr;

	auto beg = builders.begin();
	auto end = builders.end();
	while (beg != end) {
		delete (beg++)->second;
	}
	builders.clear();
}
void LocalMSGBuilder::build(const ScoreVector & score_vector) {
	/* get mutation block of the mutant */
	Mutant::ID mid = score_vector.get_mutant();
	if (!(graph->has_block_of(mid))) return;
	MutBlock & block = graph->get_block_of(mid);

	/* get the builder for local MSG in the block */
	MSGBuilder * builder;
	if (builders.count(&block) == 0) {
		builder = new MSGBuilder();
		builders[&block] = builder;
		builder->open(block.get_graph_for_modify());
	}
	else {
		auto iter = builders.find(&block);
		builder = iter->second;
	}

	/* add score vector */
	builder->add_score_vector(score_vector);
}
void LocalMSGBuilder::end() {
	auto beg = builders.begin();
	auto end = builders.end();
	while (beg != end) {
		MSGBuilder * builder = (beg++)->second;
		builder->end_score_vectors();
	}
}
void LocalMSGBuilder::construct_graphs(MutBlockGraph &bgraph, 
	ScoreProducer & producer, ScoreConsumer & consumer) {
	ScoreVector * vector;

	this->open(bgraph);
	while ((vector = producer.produce()) != nullptr) {
		this->build(*vector);
	}
	this->end(); this->close();
}

void MutBlockConnect::open(MutBlockGraph & bgraph) {
	close(); graph = &bgraph; 
	graph->clear_bridges();
}
void MutBlockConnect::close() { graph = nullptr; }
void MutBlockConnect::build_all_bridges() {
	if (graph == nullptr) return;
	const std::set<MutBlock *> blocks = graph->get_blocks();

	auto src_beg = blocks.begin(), src_end = blocks.end();
	while (src_beg != src_end) {
		/* get next source block */
		MutBlock & src = *(*(src_beg++));
		const BitSeq & src_coverage = src.get_coverage();
		if (src_coverage.all_zeros()) continue;

		auto trg_beg = blocks.begin();
		auto trg_end = blocks.end();
		while (trg_beg != trg_end) {
			/* get next unequal-nodes in graph */
			MutBlock & trg = *(*(trg_beg++));
			if ((&src) == (&trg)) continue;

			/* validate whether bridge between blocks */
			BitSeq inter = trg.get_coverage();
			inter.conjunct(src_coverage);
			if (inter.all_zeros()) continue;

			/* create new (empty) bridge between them */
			graph->new_bridge(src, trg);
		}
	}
}
void MutBlockConnect::build_all_subsume() {
	if (graph == nullptr) return;
	const std::set<MutBlockBridge *> & bridges = graph->get_bridges();
	auto beg = bridges.begin(), end = bridges.end();
	while (beg != end) {
		MutBlockBridge * bridge = *(beg++);
		this->build_subsume(*bridge);
	}
}
void MutBlockConnect::connect(MutBlockGraph & bgraph) {
	this->open(bgraph);
	this->build_all_bridges();
	this->build_all_subsume();
	this->close();
}

void MutBridgeLinker::connect(MutBlockBridge & brg) {
	this->open(brg);
	while (!vex_queue.empty()) {
		MSGVertex * vex = vex_queue.front();
		vex_queue.pop(); this->link(vex);
	}
	this->close();
}
void MutBridgeLinker::open(MutBlockBridge & brg) {
	/* cache for constructing sort-queue */
	std::queue<MSGVertex *> tqueue;
	std::set<MSGVertex *> records;

	/* initialization */
	close(); bridge = &brg;

	/* get leafs and initialize their solutions */
	const std::set<MSGVertex *> vertices = brg.get_vertices();
	auto beg = vertices.begin(), end = vertices.end();
	while (beg != end) {
		MSGVertex * vertex = *(beg++);
		if (is_leaf_vertex(vertex, vertices)) {
			/* initialize its solution for vertex */
			std::set<MSGVertex *> * solution
				= new std::set<MSGVertex *>();
			solutions[vertex] = solution;
			get_leaf_direct_subsuming(vertex, *solution);

			/* push to queue */ tqueue.push(vertex);
		}
	} /* end while: beg != end */

	/* construct sorted-queue for vertices */
	while (!tqueue.empty()) {
		/* get next valid node */
		MSGVertex * x = tqueue.front();
		records.insert(x); vex_queue.push(x);

		/* get input edges */
		const std::list<MuSubsume> & edges
			= x->get_in_port().get_edges();
		auto beg = edges.begin(), end = edges.end();

		/* update the tqueue */
		tqueue.pop();
		while (beg != end) {
			/* get next node subsuming the x */
			const MuSubsume & edge = *(beg++);
			MSGVertex * y = (MSGVertex *)(&(edge.get_source()));

			if (vertices.count(y) == 0) continue;				/* invalid */
			else if (records.count(y) > 0) continue;			/* duplicated */
			else if(is_solvable(y, records)) tqueue.push(y);	/* insert solvable */
		} /* end while: beg != end */
	} /* end while: vex_queue */

	/* return */ return;
}
void MutBridgeLinker::close() {
	while (!vex_queue.empty()) vex_queue.pop();

	auto beg = solutions.begin(), end = solutions.end();
	while (beg != end) delete (beg++)->second;
	solutions.clear(); bridge = nullptr;
}
bool MutBridgeLinker::is_leaf_vertex(MSGVertex * vex, const std::set<MSGVertex *> & vertices) {
	/* get output edges */
	const std::list<MuSubsume> & edges 
			= vex->get_ou_port().get_edges();
	auto beg = edges.begin(), end = edges.end();
	while (beg != end) {
		const MuSubsume & edge = *(beg++);
		MSGVertex * trg = (MSGVertex *)(&(edge.get_target()));
		if (vertices.count(trg) != 0) return false;
	}
	return true;
}
bool MutBridgeLinker::is_solvable(MSGVertex * x, const std::set<MSGVertex *> & records) {
	/* get output edges */
	const std::list<MuSubsume> & edges
		= x->get_ou_port().get_edges();
	auto beg = edges.begin(), end = edges.end();
	while (beg != end) {
		const MuSubsume & edge = *(beg++);
		MSGVertex * trg = (MSGVertex *)(&(edge.get_target()));
		if (records.count(trg) == 0) return false;
	}
	return true;
}
void MutBridgeLinker::get_leaf_direct_subsuming(MSGVertex * leaf, std::set<MSGVertex *> & DS) {
	/* initialization */
	const std::set<MSGVertex *> & leafs = bridge->
		get_target_block().get_local_graph().get_leafs();
	auto beg = leafs.begin(), end = leafs.end(); DS.clear();

	/* put leafs in DS */
	while (beg != end) DS.insert(*(beg++));
}
void MutBridgeLinker::link(MSGVertex * vex) {
	/* get the solution of vex for computed */
	std::set<MSGVertex *> * DS;
	if (solutions.count(vex) == 0) {
		DS = new std::set<MSGVertex *>();
		solutions[vex] = DS;
	}
	else {
		auto iter = solutions.find(vex);
		DS = iter->second;
	}

	/* compute and connect DS */
	initial_direct_subsuming(*vex, *DS);
	compute_direct_subsuming(*vex, *DS);
	connect_direct_subsuming(*vex, *DS);
}
void MutBridgeLinker::initial_direct_subsuming(MSGVertex & x, std::set<MSGVertex *> & DS) {
	/* get the output edges for x */
	const std::list<MuSubsume> & edges 
		= x.get_ou_port().get_edges();
	auto beg = edges.begin(), end = edges.end(); 
	
	/* iterate all its valid direct children */
	const std::set<MSGVertex *> & vertices = bridge->get_vertices();
	while (beg != end) {
		/* get the next target vertex of x */
		const MuSubsume & edge = *(beg++);
		MSGVertex * y = (MSGVertex *)(&(edge.get_target()));
		if (vertices.count(y) == 0) continue;

		/* get the solution DS for y */
		auto iter = solutions.find(y);
		std::set<MSGVertex *> & y_DS = *(iter->second);

		/* add y.DS to x.DS */
		auto ybeg = y_DS.begin(), yend = y_DS.end();
		while (ybeg != yend) DS.insert(*(ybeg++));
	} /* end while: beg != end */

	/* return */ return;
}
void MutBridgeLinker::compute_direct_subsuming(MSGVertex & x, std::set<MSGVertex *> & DS) {
	/* declarations */
	std::queue<MSGVertex *> vqueue;
	std::set<MSGVertex *> visits;

	/* initialize the vqueue and visit-set */
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MSGVertex & y = *(*(beg++));
		if (subsume(x, y)) {
			vqueue.push(&y);
			visits.insert(&y);
		}
	}
	DS.clear();

	/* iterate by BFS */
	while (!vqueue.empty()) {
		/* get next unvisited y in queue */
		MSGVertex & y = *(vqueue.front()); 
		vqueue.pop();
		if (visits.count(&y) > 0) continue;

		bool is_directly_subsumed = true; /* assumed y as DS */

		/* iterate y's subsuming mutations */
		const std::list<MuSubsume> & edges = y.get_in_port().get_edges();
		auto edge_beg = edges.begin(), edge_end = edges.end();
		while (edge_beg != edge_end) {
			/* get next parent (non-compared) */
			const MuSubsume & edge = *(edge_beg++);
			MSGVertex * parent = (MSGVertex *)(&(edge.get_source()));
			if (visits.count(parent) > 0) continue;

			/* insert the parent to the queue by temperary */
			if (subsume(x, *parent)) {
				is_directly_subsumed = false;
				vqueue.push(parent);
				visits.insert(parent);
			}
		} /* end while: edges iteration */

		if (is_directly_subsumed) DS.insert(&y);
	} /* end while: vqueue */

	/* return */ return;
}
void MutBridgeLinker::connect_direct_subsuming(MSGVertex & x, std::set<MSGVertex *> & DS) {
	auto beg = DS.begin(), end = DS.end();
	while (beg != end) {
		MSGVertex & y = *(*(beg++));
		bridge->links(x, y);
	}
}
bool MutBridgeLinker::subsume(MSGVertex & x, MSGVertex & y) {
	const BitSeq & xvec = x.get_feature()->get_vector();
	const BitSeq & yvec = y.get_feature()->get_vector();
	return xvec.subsume(yvec);
}

int main() {
	// initialization
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string prname = "mid";
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
		graph_builder.construct_graphs(blocks, sproducer, sconsumer);
		block_connect.connect(blocks);

		// print block information 
		std::cout << "\t[B]: " << blocks.get_blocks().size() << std::endl;
		auto beg = blocks.get_blocks().begin();
		auto end = blocks.get_blocks().end();
		while (beg != end) {
			MutBlock & block = *(*(beg++));
			std::cout << "\t\t[B]: { M = " << block.get_mutants().size();
			const MSGraph & local_graph = block.get_local_graph();
			std::cout << "; C = " << local_graph.number_of_vertices();
			std::cout << " }\n";
		}

		

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