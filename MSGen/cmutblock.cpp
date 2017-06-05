#include "cmutblock.h"

CMutantBlock * CMutantBlockSet::new_block(const BitSeq & covvec) {
	CMutantBlock * block = new CMutantBlock(covvec, mspace);
	blocks.insert(block); return block;
}
CMutantBlock & CMutantBlockSet::get_block(Mutant::ID mid) {
	if (index.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, "CMutantBlockSet::get_block", "Invalid argument mid (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error);
	}
	
	auto iter = index.find(mid);
	auto value = iter->second;
	return *value;
}
void CMutantBlockSet::add_mutant(Mutant::ID mid, CMutantBlock * block) {
	if (block == nullptr) {
		CError error(CErrorType::InvalidArguments, "CMutantBlockSet::add_mutant", "Invalid argument block (nullptr)");
		CErrorConsumer::consume(error);
	}
	else if (index.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments, "CMutantBlockSet::add_mutant", "Duplicated argument mid (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error);
	}
	else index[mid] = block;
}
void CMutantBlockSet::clear() {
	auto beg = blocks.begin();
	auto end = blocks.end();
	while (beg != end) {
		delete *(beg++);
	}
	blocks.clear();
}

void CMutantBlockBuilder::open(CMutantBlockSet & set) {
	close();
	bset = &set;
	trie = new BitTrieTree();
}
void CMutantBlockBuilder::close() {
	delete trie; 
	trie = nullptr;
	bset = nullptr;
}
void CMutantBlockBuilder::insert(const CoverageVector & covvec) {
	/* insert coverage into the trie */
	BitTrie * leaf = trie->insert_vector(covvec.get_coverage());

	/* get the block by the coverage */
	CMutantBlock * block = nullptr;
	if (leaf->get_data() == nullptr) {
		block = bset->new_block(covvec.get_coverage());
		leaf->set_data(block);
	}
	else block = (CMutantBlock *)(leaf->get_data());

	/* put mutant into the block */
	block->add_mutant(covvec.get_mutant());
	bset->add_mutant(covvec.get_mutant(), block);
}
void CMutantBlockBuilder::build_for(CMutantBlockSet & bset, 
	CoverageProducer & producer, CoverageConsumer & consumer) {
	CoverageVector * covvec;

	this->open(bset);
	while ((covvec = producer.produce()) != nullptr) {
		this->insert(*covvec); consumer.consume(covvec);
	}
	this->close();
}

void LocalMSGraphBuilder::open(CMutantBlockSet & set) {
	close(); bset = &set; 
	auto beg = set.get_blocks().begin();
	auto end = set.get_blocks().end();
	while (beg != end) {
		CMutantBlock * block = *(beg++);
		MSGBuilder * gbuilder = new MSGBuilder();
		gbuilder->open(block->get_graph_for_build());
		builders[block] = gbuilder;
	}
}
void LocalMSGraphBuilder::close() {
	/* clear builders */
	auto beg = builders.begin();
	auto end = builders.end();
	while (beg != end) {
		MSGBuilder * builder = (beg++)->second;
		builder->end_score_vectors();
		delete builder;
	}
	builders.clear();

	/* reset blocks */ bset = nullptr;
}
void LocalMSGraphBuilder::insert(const ScoreVector & svec) {
	CMutantBlock & block = bset->get_block(svec.get_mutant());
	MSGBuilder * builder = (builders.find(&block))->second;
	builder->add_score_vector(svec);
}
void LocalMSGraphBuilder::build_local_graph(CMutantBlockSet & set, 
	ScoreProducer & producer, ScoreConsumer & consumer) {
	ScoreVector * svec;

	this->open(set);
	while ((svec = producer.produce()) != nullptr) {
		this->insert(*svec);
		consumer.consume(svec);
	}
	this->close();

	return;
}

int main() {
	// initialization
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string prname = "tcas"; 
	File & root = *(new File(prefix + prname)); TestType ttype = TestType::tcas;

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
	CMutantBlockBuilder builder; LocalMSGraphBuilder gbuilder;

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
		CMutantBlockSet & blocks = *(new CMutantBlockSet(mspace));

		// building blocks
		builder.build_for(blocks, cproducer, ccomsumer);
		gbuilder.build_local_graph(blocks, sproducer, sconsumer);

		// print block information 
		std::cout << "\t[B]: " << blocks.get_blocks().size() << std::endl;
		auto beg = blocks.get_blocks().begin();
		auto end = blocks.get_blocks().end();
		while (beg != end) {
			CMutantBlock & block = *(*(beg++));
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