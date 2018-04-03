// include-list
#include "sgraph.h"
#include "cfunc.h"
#include <time.h>

/* ------------------ Basic Methods ------------------------- */
/* eliminate the spaces of '\n' */
static void trim_string(std::string & str) {
	std::string cache;
	int n = str.length();
	for (int i = 0; i < n; i++) {
		if (str[i] != '\n')
			cache += str[i];
	}
	str = cache;
}
/* ------------------ Basic Methods ------------------------- */

/* ------------------ Initialization Methods ------------------------- */
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
/* classify by functions */
static void classify_by_functions(MutantSpace & space, const CFunctionSpace & funcs, 
	std::map<Mutant::ID, std::string> & typelib) {
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
/* ------------------ Initialization Methods ------------------------- */


/* ------------------ Mutation Selection ------------------------- */
/* manually establish the operators used in experiment */
static void select_operators(std::set<std::string> & operators) {
	operators.clear();

	// incremental errors
	operators.insert("u-VDTR");
	operators.insert("u-VTWD");
	operators.insert("I-DirVarIncDec");
	operators.insert("II-ArgIncDec");
	operators.insert("I-IndVarIncDec");
	// constant errors
	operators.insert("u-Cccr");
	operators.insert("u-Ccsr");
	operators.insert("u-CRCR");
	// negative
	operators.insert("u-OBNG");
	operators.insert("u-OCNG");
	operators.insert("u-OLNG");
	operators.insert("I-DirVarAriNeg");
	operators.insert("I-DirVarBitNeg");
	operators.insert("I-DirVarLogNeg");
	operators.insert("I-IndVarAriNeg");
	operators.insert("I-IndVarBitNeg");
	operators.insert("I-IndVarLogNeg");
	// arithmetic
	operators.insert("u-OAAN");
	operators.insert("u-OABN");
	operators.insert("u-OALN");
	operators.insert("u-OARN");
	operators.insert("u-OASN");
	// bitwise
	operators.insert("u-OBAN");
	operators.insert("u-OBBN");
	operators.insert("u-OBLN");
	operators.insert("u-OBRN");
	operators.insert("u-OBSN");
	// logical connector
	operators.insert("u-OLAN");
	operators.insert("u-OLBN");
	operators.insert("u-OLLN");
	operators.insert("u-OLRN");
	operators.insert("u-OLSN");
	// relational
	operators.insert("u-ORAN");
	operators.insert("u-ORBN");
	operators.insert("u-ORLN");
	operators.insert("u-ORRN");
	operators.insert("u-ORSN");
	// shifting error
	operators.insert("u-OSAN");
	operators.insert("u-OSBN");
	operators.insert("u-OSLN");
	operators.insert("u-OSRN");
	operators.insert("u-OSSN");
	// statement
	operators.insert("u-SMTC");
	operators.insert("u-SSDL");
	operators.insert("u-STRI");
	operators.insert("u-STRP");
	// variable
	operators.insert("u-VLSR");
	operators.insert("u-VGSR");

}
/* select mutations by operators */
static void select_mutants(MutantSpace & mspace, 
	const std::set<std::string> & operators, 
	std::set<Mutant::ID> & mutants) {

	mutants.clear();	// clear all remainders

	// select mutant by operator
	Mutant::ID n = mspace.number_of_mutants();
	for (Mutant::ID mid = 0; mid < n; mid++) {
		Mutant & mk = mspace.get_mutant(mid);
		if (operators.count(mk.get_operator()) > 0) {
			mutants.insert(mid);
		}
	}

	return;		// ending
}
/* ------------------ Mutation Selection ------------------------- */


/* ------------------ MSG Building ------------------------- */
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
/* build up MSG from mutants in specified subset of mutants */
static void load_ms_graph(
	const File & root, const CodeFile & cfile, 
	CMutant & cmutant, CTest & ctest, CScore & cscore, 
	const std::set<Mutant::ID> & subset, 
	MS_Graph & graph) {
	// get set of mutants and tests in project
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();

	// get score vector producer | consumer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	FileScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	// filter for selection
	ScoreFilter fproducer(producer, subset);
	build_up_graph(graph, fproducer, consumer);

	// release resource
	ctest.delete_test_set(&tests); mspace.delete_set(&mutants);
}
/* build up MSG from project data */
static void load_ms_graph(
	const File & root, const CodeFile & cfile, 
	CMutant & cmutant, CTest & ctest, CScore & cscore, 
	MS_Graph & graph) {
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
/* ------------------ MSG Building ------------------------- */


/* ------------------ Analysis Methods ------------------------- */
/* get the nodes subsumed by specified node */
static void get_subsumed(
	MSG_Node & root, 
	std::set<MSG_Node *> & children) {
	// initialization
	children.clear();
	std::queue<MSG_Node *> queue;
	queue.push(&root);

	// iterate the children 
	MSG_Node * node;
	while (!queue.empty()) {
		node = queue.front();
		queue.pop();

		const MSG_Port & port 
			= node->get_ou_port();
		int n = port.degree();

		for (int i = 0; i < n; i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();

			if (children.count(&next) == 0) {
				children.insert(&next);
				queue.push(&next);
			}
		} // end for

	} // end while

	return;		// end
}
/* get the nodes that subsume the leaf */
static void get_subsuming(
	MSG_Node & leaf,
	std::set<MSG_Node *> & parents) {
	// initialization
	parents.clear();
	std::queue<MSG_Node *> queue;
	queue.push(&leaf);

	// iterate the parents
	MSG_Node * node;
	while (!queue.empty()) {
		node = queue.front();
		queue.pop();

		const MSG_Port & port 
			= node->get_in_port();
		int n = port.degree();

		for (int i = 0; i < n; i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & prev = edge.get_source();

			if (parents.count(&prev) == 0) {
				parents.insert(&prev);
				queue.push(&prev);
			}
		}
	} // end while

	return;		// end
}
/* compute the utility of a mutant node */
static double get_utility_of(MSG_Node & node) {
	std::set<MSG_Node *> nodes;
	int parents, children;

	get_subsuming(node, nodes);
	parents = nodes.size();
	get_subsumed(node, nodes);
	children = nodes.size();

	return ((double)children) 
		/ ((double)(parents + children));
}
/* ------------------ Analysis Methods ------------------------- */


/* ------------------ Outputters ------------------------- */
/* ( mid | cid | operator | function | line | original | replace ) */
static void output_mutants(
	const MS_Graph & graph,
	const std::map<Mutant::ID, std::string> & fmap,
	std::ostream & out) {
	// declarations
	MutantSpace & mspace = graph.get_space();
	const TextBuild * text = mspace.get_code_file().get_text();

	// title print
	out << "mid\tcid\toperator\tfunction\tline\toriginal\treplace\n";

	// iterate 
	Mutant::ID n = mspace.number_of_mutants();
	for (Mutant::ID mid = 0; mid < n; mid++) {
		// filter
		if (!graph.has_node_of(mid)) 
			continue;
		// print mutant information
		else {
			// get information in mutant
			Mutant & mutant = mspace.get_mutant(mid);
			MSG_Node & node = graph.get_node_of(mid);
			const Mutation & mutation = mutant.
				get_mutation(mutant.get_orders() - 1);
			const CodeLocation & loc = mutation.get_location();
			std::string origin = mutation.get_location().get_text_at();
			std::string replace = mutation.get_replacement();
			trim_string(origin); trim_string(replace);
			
			// output 
			out << mid << "\t";
			out << node.get_node_id() << "\t";
			out << mutant.get_operator() << "\t";
			if (fmap.count(mid) == 0)
				out << "??????" << "\t";
			else {
				auto iter = fmap.find(mid);
				out << iter->second << "\t";
			}
			out << text->lineOfIndex(loc.get_bias()) << "\t";
			out << origin << "\t";
			out << replace << "\t";

			// print line
			out << "\n";
		}
	} // end for

	// return
	out << std::endl; return;
}
/* node | degree | mutants | nextset */
static void output_graphic(const MS_Graph & graph, std::ostream & out) {
	// title
	out << "id\tmutants\tdegree\tutility\tnextset\tdnodes\tdmutants\tcnodes\tcmutants\n";

	long n = graph.size(), k;
	std::set<MSG_Node *> nodes;
	for (k = 0; k < n; k++) {
		// get next mutant
		MSG_Node & node = graph.get_node(k);
		if (node.get_score_degree() == 0)
			continue;

		// basic information
		out << node.get_node_id() << "\t";
		out << node.get_mutants().number_of_mutants() << "\t";
		out << node.get_score_degree() << "\t";
		out << get_utility_of(node) << "\t";

		// print next nodes
		int dnodes = 0, dmutants = 0;
		const MSG_Port & port = node.get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();

			out << next.get_node_id() << "; ";

			dnodes += 1;
			dmutants += next.get_mutants().number_of_mutants();
		}
		out << "\t";

		// direct subsumption
		out << dnodes << "\t" << dmutants << "\t";

		// output other subsumption
		get_subsumed(node, nodes);
		int ch_num = 0;
		auto beg = nodes.begin();
		auto end = nodes.end();
		while (beg != end) 
			ch_num += (*(beg++))->get_mutants().number_of_mutants();
		out << nodes.size() << "\t";
		out << ch_num << "\t";

		// output line
		out << '\n';
	}
	out << std::endl;
}
/* ------------------ Outputters ------------------------- */

/* ------------------ Printters ------------------------- */
static void print_msg_prevalence(const MS_Graph & graph) {
	/* declarations */
	int n = graph.size(), m = 0, e = 0;
	unsigned long all_equiv = 0, ess_equiv = 0;
	unsigned long all_stric = 0, ess_stric = 0;
	std::set<MSG_Node *> children;

	/* compute the number of relationships */
	for (int i = 0; i < n; i++) {
		MSG_Node & node = graph.get_node(i);
		if (node.get_score_degree() > 0) {
			int me = node.get_mutants().number_of_mutants();
			int se = node.get_ou_port().degree();
			get_subsumed(node, children);

			m = m + me; e = e + se;
			ess_equiv = ess_equiv + (me - 1);
			all_equiv = all_equiv + (me - 1) * me;
			ess_stric = ess_stric + se;
			all_stric = all_stric + children.size() - 1;
		}
	}
	std::cout << "\n================ MSG-Prevalence =============\n";
	std::cout << "\tMutants\tNodes\tEdges\tAll-Equiv\tEss-Equiv\tAll-Strict\tEss-Strict\n";
	std::cout << "\t" << m << "\t" << n << "\t" << e << "\t" 
		<< all_equiv << "\t" << ess_equiv << "\t" << all_stric << "\t" << ess_stric << "\n";
	std::cout << "================ MSG-Prevalence =============\n";
}
/* ------------------ Printters ------------------------- */

/* ------------------ Main Tester ------------------------- */
int main(int argc, char *argv[]) {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "triangle";
	std::string pdir = prefix + prname + "/analysis/";
	TestType ttype = TestType::general;

	// get root file and analysis dir 
	File & root = *(new File(prefix + prname));

	// create code-project, mutant-project, test-project
	CProgram & program = *(new CProgram(root));
	CFuncProject & funcs = *(new CFuncProject(root));
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
		// get code file
		const CodeFile & cfile = *(*(beg++));
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);

		// load functions
		std::map<Mutant::ID, std::string> fmap;
		load_functions(funcs, cmutant, cfile, fmap);

		// test procedure
		MS_Graph graph(mspace);
		load_ms_graph(root, cfile, cmutant, ctest, cscore, graph);

		// output procedure
		std::ofstream out1(pdir + "mutants.txt");
		output_mutants(graph, fmap, out1); out1.close();
		std::ofstream out2(pdir + "msgraph.txt");
		output_graphic(graph, out2); out2.close();

		// print-procedure
		print_msg_prevalence(graph);

		// end this file
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &cscore; delete &funcs;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	// end of all
	std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
