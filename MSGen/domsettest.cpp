#include "domset.h"
#include "cfunc.h"
#include <time.h>

// tag methods
/* classify by operator */
static void classify_by_operators(MutantSpace & space, std::map<Mutant::ID, std::string> & typelib) {
	typelib.clear();
	Mutant::ID mid, n = space.number_of_mutants();
	for (mid = 0; mid < n; mid++) {
		Mutant & mutant = space.get_mutant(mid);
		typelib[mutant.get_id()] = mutant.get_operator();
	}
}
/* classify by functions */
static void classify_by_functions(MutantSpace & space, const CFunctionSpace & funcs, std::map<Mutant::ID, std::string> & typelib) {
	typelib.clear(); std::string funcname;
	Mutant::ID mid, n = space.number_of_mutants();
	for (mid = 0; mid < n; mid++) {
		/* get the mutant's location */
		Mutant & mutant = space.get_mutant(mid);
		const Mutation & mutation = mutant.
			get_mutation(mutant.get_orders() - 1);
		const CodeLocation & loc = mutation.get_location();

		if (funcs.find_function_at(loc.get_bias(), funcname))
			typelib[mid] = funcname;
		//else typelib[mid] = "unknown";
	}
}

// basic methods
/* load the tests and mutants into the project */
static void load_tests_mutants(CTest & ctest, CMutant & cmutant) {
	ctest.load(); const TestSpace & tspace = ctest.get_space();
	std::cout << "Loading test cases: " << tspace.number_of_tests() << std::endl;

	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		/* load file text */
		CodeFile & cfile = *(*(cfile_beg++));
		cfile.get_space().load(cfile);
		/* load mutants and mutations */
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		cmutant.load_mutants_for(mspace, true);
		std::cout << "Load " << mspace.number_of_mutants() <<
			" mutants for: " << cfile.get_file().get_path() << "\n" << std::endl;
	}
}
/* load the functions */
static void load_functions(CFuncProject & funcs, CMutant & cmutant,
	const CodeFile & cfile, std::map<Mutant::ID, std::string> & typelib) {
	funcs.load_functions_for(cfile);
	CFunctionSpace & funcspace = funcs.get_function_space(cfile);

	// load functions
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	classify_by_functions(mspace, funcspace, typelib);

	// head information
	std::cout << "Load functions:\n\t{";
	const std::set<std::string> & funcnames = funcspace.get_function_names();
	auto beg = funcnames.begin(), end = funcnames.end();
	while (beg != end) {
		std::cout << *(beg++);
		if (beg != end) std::cout << "; ";
	}
	std::cout << "}\n";
}
/* load score set to score matrix */
static ScoreMatrix & load_score_set(const CodeFile & cfile, 
	CMutant & cmutant, CTest & ctest, CScore & cscore) {
	// get set of mutants and tests in project
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();

	// get score vector producer | consumer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	FileScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	ScoreMatrix * matrix = new ScoreMatrix(mspace, ctest.get_space());
	matrix->add_score_vectors(producer, consumer);

	// release resource
	ctest.delete_test_set(&tests); mspace.delete_set(&mutants);

	return *matrix;		// return dynamic matrix
}

// file-output methods
static void output_testing(const std::string & file, DomSetBuilder & builder) {
	// declarations
	ScoreMatrix & input = builder.get_input_matrix();
	MutantSpace & mspace = input.get_mutant_space();
	Mutant::ID n = mspace.number_of_mutants();
	DomSetAlgorithm_Counter & counter
		= builder.get_performance_counter();

	// outputs
	std::ofstream out(file);
	out << "Mutant\t#test\n";
	for (Mutant::ID mid = 0; mid < n; mid++) {
		size_t tnum = counter.get_mutant_tests(mid);
		out << mid << "\t" << tnum << "\n";
	}
	out << std::endl; out.close();
}
static void output_compare(const std::string & file, DomSetBuilder & builder) {
	// declarations
	ScoreMatrix & input = builder.get_input_matrix();
	MutantSpace & mspace = input.get_mutant_space();
	Mutant::ID n = mspace.number_of_mutants();
	DomSetAlgorithm_Counter & counter
		= builder.get_performance_counter();

	std::ofstream out(file);
	out << "Mutant\t#trans\n";
	for (Mutant::ID mid = 0; mid < n; mid++) {
		size_t tnum = counter.get_states_trans(mid);
		out << mid << "\t" << tnum << "\n";
	}
	out << std::endl; out.close();
}
static void output_decline(const std::string & file, DomSetBuilder & builder) {
	// declarations
	DomSetAlgorithm_Counter & counter
		= builder.get_performance_counter();

	std::ofstream out(file);
	out << "Mutant\t#DSet\n";

	const std::vector<MutPair> &
		declines = counter.get_declines();
	for (int i = 0; i < declines.size(); i++) {
		MutPair pair = declines[i];
		out << pair.mid << "\t" << pair.declines << "\n";
	}

	out << std::endl; out.close();
}

// compute methods
/* build up MSG */
static void build_up_graph(MS_Graph & graph, ScoreProducer & producer, ScoreConsumer & consumer) {
	time_t start, end;
	//MSG_Build & builder = *(new MSG_Build_Classical(graph));
	MSG_Build & builder = *(new MSG_Build_Fast(graph));
	//MSG_Build & builder = *(new MSG_Build_Quick(graph));
	builder.open(producer, consumer);
	start = time(nullptr);
	builder.build();
	end = time(nullptr);
	builder.close();
	delete &builder;
	std::cout << "\tUsing time: " << difftime(end, start) << " seconds...\n";
}
/* build up MSG from project data */
static void load_ms_graph(const File & root, const CodeFile & cfile, CFuncProject & funcs,
	CMutant & cmutant, CTest & ctest, CScore & cscore, MS_Graph & graph) {
	// get set of mutants and tests in project
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();

	// get score vector producer | consumer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	FileScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	// MS-Graph-Build
	build_up_graph(graph, producer, consumer);

	// release resource
	ctest.delete_test_set(&tests); mspace.delete_set(&mutants);
}
/* build up CSG from project data */
static void load_cs_graph(const File & root, const CodeFile & cfile,
	CMutant & cmutant, CTest & ctest, CTrace & ctrace, CScore & cscore, MS_Graph & graph) {
	// load coverage information
	CoverageSpace & covspac = ctrace.get_space();
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
	covspac.add_file_coverage(cfile, tests); ctrace.load_coverage(cfile);
	std::cout << "Coverage loading...\n";

	// coverage vector producer
	FileCoverage & filecov = covspac.get_file_coverage(cfile);
	CoverageProducer cproducer(mspace, filecov);
	CoverageConsumer cconsumer;

	// score producer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	CoverageScoreProducer producer(score_func, cproducer, cconsumer);
	ScoreConsumer consumer(score_func);

	// MS-Graph-Build
	build_up_graph(graph, producer, consumer);

	// release resource
	ctest.delete_test_set(&tests);
}
/* compute dominator set by greedy algorithm */
static void compute_dominator_set_by_greedy(ScoreMatrix & matrix, MutSet & domset, File & root) {
	// compute the dominator set by greedy algorithm
	DomSetBuilder_Greedy builder;
	time_t start = time(nullptr);
	builder.open(matrix);
	builder.build(domset);
	time_t end = time(nullptr);

	// efficiency analysis
	std::string prefix = root.get_path() + "/analysis/";
	output_testing(prefix + "greedy_testing.txt", builder);
	output_compare(prefix + "greedy_compare.txt", builder);
	output_decline(prefix + "greedy_decline.txt", builder);

	// end
	builder.close();
	std::cout << "/------ Greedy Algorithm -------/\n";
	std::cout << "\tEquivalent: " << matrix.get_equivalents() << "\n";
	std::cout << "\tSubsumings: " << domset.size() << "\n";
	std::cout << "/-------------------------------/\n";
}
/* compute dominator set by coverage algorithm */
static void compute_dominator_set_by_coverage(ScoreMatrix & matrix, MS_Graph & csg, MutSet & domset, File & root) {
	// compute dominator set by coverage algorithm
	DomSetBuilder_Blocks builder(csg);
	time_t start = time(nullptr);
	builder.open(matrix);
	builder.build(domset);
	time_t end = time(nullptr);

	// efficiency analysis
	std::string prefix = root.get_path() + "/analysis/";
	output_testing(prefix + "covers_testing.txt", builder);
	output_compare(prefix + "covers_compare.txt", builder);
	output_decline(prefix + "covers_decline.txt", builder);

	// end
	builder.close();
	std::cout << "/------ Optimized Algorithm -------/\n";
	std::cout << "\tEquivalent: " << matrix.get_equivalents() << "\n";
	std::cout << "\tSubsumings: " << domset.size() << "\n";
	std::cout << "/-------------------------------/\n";
}

// main method
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "replace";
	TestType ttype = TestType::replace;

	// get root file and analysis dir 
	File & root = *(new File(prefix + prname));

	// create code-project, mutant-project, test-project
	CProgram & program = *(new CProgram(root));
	CFuncProject & funcs = *(new CFuncProject(root));
	CTest & ctest = *(new CTest(ttype, root, program.get_exec()));
	CMutant & cmutant = *(new CMutant(root, program.get_source()));
	CScore & cscore = *(new CScore(root, cmutant, ctest));
	CTrace & ctrace = *(new CTrace(root, program.get_source(), ctest.get_space()));

	// load mutants and tests
	load_tests_mutants(ctest, cmutant);

	// dominator set for each file 
	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto beg = cfiles.begin(), end = cfiles.end();
	while (beg != end) {
		// get the next file of code
		const CodeFile & cfile = *(*(beg++));
		std::map<Mutant::ID, std::string> funclib;
		load_functions(funcs, cmutant, cfile, funclib);
		std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

		// load csg 
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		MS_Graph cgraph(mspace);
		load_cs_graph(root, cfile, cmutant, ctest, ctrace, cscore, cgraph);

		// get dominator set
		MutSet domset(mspace); 
		ScoreMatrix & matrix = load_score_set(cfile, cmutant, ctest, cscore);
		std::cout << "Equivalent mutants: " << matrix.get_equivalents() << "\n";
		compute_dominator_set_by_greedy(matrix, domset, root);
		compute_dominator_set_by_coverage(matrix, cgraph, domset, root);
		delete &matrix;

		// end this file
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	} // end while

	// delete memory
	delete &ctrace;
	delete &cscore; delete &funcs;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
