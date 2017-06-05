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

int main() {
	// initialization
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string prname = "profit"; 
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

	// get coverage vector
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
	cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		// get next file and load its text 
		CodeFile & cfile = *(*(cfile_beg++)); cspace.load(cfile);

		// get mutations for code-file
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);

		// load coverage 
		covspace.add_file_coverage(cfile, tests);
		ctrace.load_coverage(cfile);
		std::cout << "Loading coverage for \"" << cfile.get_file().get_path() << "\"\n";

		// get coverage producer and consumer
		FileCoverage & fcov = covspace.get_file_coverage(cfile);
		CoverageProducer cproducer(mspace, fcov);
		CoverageConsumer ccomsumer;

		// create mutblock and builder
		CMutantBlockSet & blocks = *(new CMutantBlockSet(mspace));
		CMutantBlockBuilder builder;		

		// building blocks
		builder.build_for(blocks, cproducer, ccomsumer);
		std::cout << "\t[B]: " << blocks.get_blocks().size() << std::endl;
		auto beg = blocks.get_blocks().begin();
		auto end = blocks.get_blocks().end();
		while (beg != end) {
			CMutantBlock & block = *(*(beg++));
			std::cout << "\t\t[C]:" << block.get_mutants().size() << "\t" << block.get_coverage().to_string() << std::endl;
		}

		// delete resources
		delete &blocks;
	}

	// delete resources
	delete &ctrace;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;
	// end to return 
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}