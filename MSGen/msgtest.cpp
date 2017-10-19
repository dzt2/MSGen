#include "sgraph.h"
#include <time.h>

typedef std::string MClass;

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
/* build up MSG */
static void build_up_graph(MS_Graph & graph, ScoreProducer & producer, ScoreConsumer & consumer) {
	time_t start, end;
	//MSG_Build & builder = *(new MSG_Build_Classical(graph));
	//MSG_Build & builder = *(new MSG_Build_Fast(graph));
	MSG_Build & builder = *(new MSG_Build_Quick(graph));
	builder.open(producer, consumer);
	start = time(nullptr);
	builder.build();
	end = time(nullptr);
	builder.close();
	delete &builder;
	std::cout << "\tUsing time: " << difftime(end, start) << " seconds...\n";
}
/* print the information about MSG-structure */
static void print_ms_graph(const MS_Graph & graph, std::ostream & out) {
	// output structure summary
	int C = graph.size(), S = 0, E = 0;
	int M = graph.size_of_mutants();
	for (int k = 0; k < C; k++) {
		MSG_Node & node = graph.get_node(k);
		S += node.get_ou_port().degree();
		if (node.get_score_degree() == 0)
			E = node.get_mutants().number_of_mutants();
	}
	int R = (M - E) * (M - E) - M;
	out << "---------- Summary ----------\n";
	out << "\tMutants: " << M << "\n";
	out << "\tEquivnt: " << E << "\n";
	out << "\tRelations: " << R << "\n";
	out << "-----------------------------\n\n";

	// total relations
	out << "---------- Graph ----------\n";
	out << "\tClusters: " << C << "\n";
	out << "\tSubsumes: " << S << "\n";

	// output equivalence edges
	int ER = 0, Ci, EDR = 0;
	for (int i = 0; i < C; i++) {
		MSG_Node & node = graph.get_node(i);
		if (node.get_score_degree() > 0) {
			Ci = node.get_mutants().number_of_mutants();
			ER += Ci * (Ci - 1) / 2;
		}
		else EDR = node.get_ou_port().degree();
	}
	double erate = 2.0 * ((double)ER) / ((double)R);
	out << "\tEqualwith: " << ER << " (" << erate << ")\n";

	// output direct subsumption
	double drate = ((double)(S - EDR)) / ((double)R);
	out << "\tDirectSub: " << (S - EDR) << " (" << drate << ")\n";
	out << "-----------------------------\n\n";

	// output eof
	out << std::endl;
}

// analysis process
/* whether the mutant is equivalent */
static bool is_equivalent_mutant(Mutant::ID i, const MS_Graph & graph) {
	if (graph.has_node_of(i)) {
		MSG_Node & node = graph.get_node_of(i);
		return node.get_score_degree() == 0;
	}
	else return false;
}
/* whether the mi and mj are indistinguishable */
static bool indistinguishable(Mutant::ID i, Mutant::ID j, const MS_Graph & graph) {
	if (graph.has_node_of(i)) {
		if (graph.has_node_of(j)) {
			MSG_Node & mi = graph.get_node_of(i);
			MSG_Node & mj = graph.get_node_of(j);
			return (&mi) == (&mj);
		}
		else return false;
	}
	else return false;
}

/* relation of mutants' operators */
MClass relation_of_operators(Mutant::ID i, Mutant::ID j, MutantSpace & space) {
	Mutant & mi = space.get_mutant(i);
	Mutant & mj = space.get_mutant(j);
	return "(" + mi.get_operator() 
		+ ", " + mj.get_operator() + ")";
}

int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "Day";
	TestType ttype = TestType::general;

	// get root file and analysis dir 
	File & root = *(new File(prefix + prname));

	// create code-project, mutant-project, test-project
	CProgram & program = *(new CProgram(root));
	CTest & ctest = *(new CTest(ttype, root, program.get_exec()));
	CMutant & cmutant = *(new CMutant(root, program.get_source()));
	CScore & cscore = *(new CScore(root, cmutant, ctest));

	// load mutants and tests
	load_tests_mutants(ctest, cmutant);

	// load MSG
	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto beg = cfiles.begin(), end = cfiles.end();
	while (beg != end) {
		// get set of mutants and tests in project
		const CodeFile & cfile = *(*(beg++));
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		MutantSet & mutants = *(mspace.create_set()); mutants.complement();
		TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
		std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

		// get score vector producer | consumer
		ScoreSource & score_src = cscore.get_source(cfile);
		ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
		ScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

		// MS-Graph-Build
		MS_Graph graph(mspace);
		build_up_graph(graph, producer, consumer);

		// print graph
		print_ms_graph(graph, std::cout);

		// end this file
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &cscore;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
