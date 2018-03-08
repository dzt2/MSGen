#include "sgraph.h"
#include "cfunc.h"
#include <time.h>

typedef std::string MClass;
static const double category_1_limie = 0.90;

// basic methods
/* trim the line-space */
static void trim_spaces(std::string & line) {
	std::string cache; char ch;
	for (int i = 0; i < line.length(); i++) {
		ch = line[i];
		if (ch != '\n') cache += ch;
	}
	line = cache;
}
/* category-method */
static bool is_category_1(Mutant::ID mid, 
	const MS_Graph & msg, const MS_Graph & csg) {
	if (msg.has_node_of(mid) && csg.has_node_of(mid)) {
		MSG_Node & mnode = msg.get_node_of(mid);
		MSG_Node & cnode = csg.get_node_of(mid);
		const BitSeq & mi = mnode.get_score_vector();
		const BitSeq & ci = cnode.get_score_vector();
		return ci.subsume(mi);
	}
	else return false;
}

// load-select data interfaces
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
/* select mutants for analysis */
static void select_mutants_by_operators(const MutantSpace & space,
	const std::set<std::string> & operators, std::set<Mutant::ID> & mutants) {
	// initialize 
	mutants.clear();

	// select mutants by operators
	Mutant::ID mid, n = space.number_of_mutants();
	for (mid = 0; mid < n; mid++) {
		Mutant & mutant = space.get_mutant(mid);
		if (operators.count(mutant.get_operator()) > 0)
			mutants.insert(mid);
	}
}
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
/* select mutants by specified function name based on its template */
static void select_mutants_by_function(
	std::map<Mutant::ID, std::string> & templt, 
	const std::string & function, 
	std::set<Mutant::ID> & mutants) {
	mutants.clear();		/* initialization */

	auto beg = templt.begin();
	auto end = templt.end();
	while (beg != end) {
		Mutant::ID mid = beg->first;
		const std::string & func = beg->second;

		if (function == func) {
			mutants.insert(mid);
		}

		beg++;
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
/* filter mutants in specified operators */
static void filter_mutants_by_operators(MutantSpace & mspace, 
	const std::set<std::string> & operators,
	const std::set<Mutant::ID> & source, 
	std::set<Mutant::ID> & target) {
	target.clear();

	auto beg = source.begin();
	auto end = source.end();
	while (beg != end) {
		Mutant::ID mid = *(beg++);
		const Mutant & mutant = mspace.get_mutant(mid);
		const std::string oprt = mutant.get_operator();
		
		if (operators.count(oprt) > 0) {
			target.insert(mutant.get_id());
		}
	}
}
/* select operators */
static void select_operators(std::set<std::string> & operators) {
	operators.clear();

	/* operator mutations */
	operators.insert("u-OAAN");
	operators.insert("u-OABN");
	operators.insert("u-OALN");
	operators.insert("u-OARN");
	operators.insert("u-OASN");

	operators.insert("u-ORAN");
	operators.insert("u-ORBN");
	operators.insert("u-ORLN");
	operators.insert("u-ORRN");
	operators.insert("u-ORSN");

	operators.insert("u-OLAN");
	operators.insert("u-OLBN");
	operators.insert("u-OLLN");
	operators.insert("u-OLRN");
	operators.insert("u-OLSN");

	operators.insert("u-OCNG");
	operators.insert("u-OLNG");

	/* key operators */
	operators.insert("u-VDTR");
	operators.insert("u-VTWD");
	//operators.insert("I-DirVarIncDec");
	//operators.insert("I-IndVarIncDec");
}

// build methods
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
/* build up MSG from mutants in specified subset of mutants */
static void load_ms_graph(const File & root, const CodeFile & cfile, CFuncProject & funcs, 
	CMutant & cmutant, CTest & ctest, CScore & cscore, const std::set<Mutant::ID> & subset, MS_Graph & graph) {
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


// calculate methods
/* get the number of nodes that subsume the leaf node */
static size_t calculate_subsuming_nodes(MSG_Node & leaf, std::set<MSG_Node *> subsumings) {
	subsumings.clear();

	std::queue<MSG_Node *> qlist;
	qlist.push(&leaf); MSG_Node * node;

	while (!qlist.empty()) {
		node = qlist.front(); qlist.pop();
		if (node->get_score_degree() > 0)
			subsumings.insert(node);

		const MSG_Port & port = node->get_in_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & prev = edge.get_source();
			if (subsumings.count(&prev) == 0) {
				qlist.push(&prev);
			}
		}
	}

	subsumings.erase(&leaf);
	return subsumings.size();
}
/* get the number of nodes that are subsumed by root */
static size_t calculate_subsummed_nodes(MSG_Node & root, std::set<MSG_Node *> & subsummeds) {
	subsummeds.clear();

	std::queue<MSG_Node *> qlist;
	qlist.push(&root); MSG_Node * node;

	while (!qlist.empty()) {
		node = qlist.front(); qlist.pop();
		if (node->get_score_degree() > 0)
			subsummeds.insert(node);

		const MSG_Port & port = node->get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();
			if (subsummeds.count(&next) == 0) {
				qlist.push(&next);
			}
		}
	}

	subsummeds.erase(&root);
	return subsummeds.size();
}
/* calculate the utility of each node in MSG */
static void calculate_mutant_utility(const MS_Graph & graph, std::map<MSG_Node *, double> & utilities) {
	utilities.clear();

	long n = graph.size(), k;
	for (k = 0; k < n; k++) {
		MSG_Node & node = graph.get_node(k);
		if (node.get_score_degree() == 0)
			utilities[&node] = 0.0;
		else {
			std::set<MSG_Node *> subsumings;
			std::set<MSG_Node *> subsummeds;
			size_t S = calculate_subsuming_nodes(node, subsumings);
			size_t D = calculate_subsummed_nodes(node, subsummeds);

			double utility; 
			if (D == 0) utility = 0.0;
			else utility = ((double)D) / ((double)(S + D));
			utilities[&node] = utility;
		}
	}
}

// print-methods
/* node | degree | mutants | nextset in dot language */
static void print_msgraph(const MS_Graph & graph, std::ostream & out) {
	/* head */ out << "digraph G {\n";

	std::map<MSG_Node *, double> utilities;
	calculate_mutant_utility(graph, utilities);

	long n = graph.size(), k;
	for (k = 0; k < n; k++) {
		MSG_Node & node = graph.get_node(k);

		/* definition of node */
		out << "n" << k << " [";
		out << "shape=box,";
		out << "label=\""; 
		out << "<MSG-Node-" << k << ">\\n";
		out << node.get_mutants().number_of_mutants() << "-mutants\\n";
		out << node.get_score_degree() << "-degrees\\n";
		auto iter = utilities.find(&node);
		out << iter->second << "-utility\\n";
		out << "\""; out << "];\n";

		/* edges information */
		const MSG_Port & port = node.get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();

			out << "n" << node.get_node_id() << " -> " 
				<< "n" << next.get_node_id() << ";\n";
		}
	}

	/* end */ out << "}" << std::endl;
}
/* id | msg_node | operator | function | line | origin | replace */
static void print_mutants(const MS_Graph & graph,
	const std::map<Mutant::ID, std::string> & funcs, std::ostream & out) {
	/* title */
	out << "id\tmsg_node\toperator\tutility\tfunction\tline\torigin\treplace\n";

	std::map<MSG_Node *, double> utilities;
	calculate_mutant_utility(graph, utilities);

	/* print mutants */
	MutantSpace & mspace = graph.get_space();
	Mutant::ID n = mspace.number_of_mutants();
	for (Mutant::ID mid = 0; mid < n; mid++) {
		/* not used mutants are ignored */
		if (!graph.has_node_of(mid)) continue;

		/* id of mutant */
		Mutant & mutant = mspace.get_mutant(mid);
		out << mutant.get_id() << "\t";

		/* msg-node id */
		MSG_Node & node = graph.get_node_of(mid);
		out << node.get_node_id() << "\t";

		/* mutant operator */
		out << mutant.get_operator() << "\t";

		/* utility */
		auto iter = utilities.find(&node);
		out << iter->second << "\t";

		/* function */
		if (funcs.count(mid) > 0) {
			auto iter = funcs.find(mid);
			out << iter->second << "\t";
		}
		else out << "unknown?\t";

		/* line */
		size_t order = mutant.get_orders();
		const Mutation & mutation = mutant.get_mutation(order - 1);
		const CodeLocation & loc = mutation.get_location();
		const TextBuild * text = loc.get_file().get_text();
		out << text->lineOfIndex(loc.get_bias()) << "\t";

		/* orign & replace */
		std::string origin = loc.get_text_at();
		std::string replace = mutation.get_replacement();
		trim_spaces(origin); trim_spaces(replace);
		out << origin << "\t" << replace << "\t";

		/* end of line */ out << "\n";
	}

	/* end */ out << std::endl;
}
/* node | degree | mutants | nextset */
static void print_graphic(const MS_Graph & graph, std::ostream & out) {
	out << "node\tdegree\tutility\tmutants\tnext-set\n";

	std::map<MSG_Node *, double> utilities;
	calculate_mutant_utility(graph, utilities);

	long n = graph.size(), k;
	for (k = 0; k < n; k++) {
		MSG_Node & node = graph.get_node(k);
		out << node.get_node_id() << "\t";
		out << node.get_score_degree() << "\t";
		auto iter = utilities.find(&node);
		out << iter->second << "\t";
		out << node.get_mutants().number_of_mutants() << "\t";

		const MSG_Port & port = node.get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();

			out << next.get_node_id() << "; ";
		}
		
		out << '\n';
	}
	out << std::endl;
}

// test module for each code file mutants
/* construct MSG for all mutants and output its information */
static void test_all_mutants_msg(const File & root, const CodeFile & cfile, 
	CFuncProject & funcs, CMutant & cmutant, CTest & ctest, CScore & cscore) {
	// get the next file of code
	std::map<Mutant::ID, std::string> funclib;
	load_functions(funcs, cmutant, cfile, funclib);
	std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

	// analysis declarations
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MS_Graph mgraph(mspace);

	// mutant subsumption graph construction 
	load_ms_graph(root, cfile, funcs, cmutant, ctest, cscore, mgraph);

	// output information
	std::ofstream out1(root.get_path() + "/analysis/all_mutants/mutants.txt");
	print_mutants(mgraph, funclib, out1); out1.close();
	std::ofstream out2(root.get_path() + "/analysis/all_mutants/msgraph.dot");
	print_msgraph(mgraph, out2); out2.close();
	std::ofstream out3(root.get_path() + "/analysis/all_mutants/msgraph.txt");
	print_graphic(mgraph, out3); out3.close();
}
/* construct MSG for specific function and output its information in matrix */
static void test_sub_mutants_msg(const File & root, const CodeFile & cfile, 
	const std::string function, CFuncProject & funcs, CMutant & cmutant, 
	CTest & ctest, CScore & cscore) {
	// get the next file of code
	std::map<Mutant::ID, std::string> funclib;
	load_functions(funcs, cmutant, cfile, funclib);
	std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

	// analysis declarations
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MS_Graph mgraph(mspace);

	// construct mutant subsumption in specified function
	std::set<Mutant::ID> source, target;
	select_mutants_by_function(funclib, function, source);
	std::set<std::string> operators; select_operators(operators);
	filter_mutants_by_operators(mspace, operators, source, target);
	std::cout << "Load mutants for \"" << function << "\": " << target.size() << "\n";
	load_ms_graph(root, cfile, funcs, cmutant, ctest, cscore, target, mgraph);

	// output graph information 
	std::ofstream out1(root.get_path() + "/analysis/all_functions/" + function + "_msgraph.dot");
	print_msgraph(mgraph, out1); out1.close();
	std::ofstream out2(root.get_path() + "/analysis/all_functions/" + function + "_mutants.txt");
	print_mutants(mgraph, funclib, out2); out2.close();
	std::ofstream out3(root.get_path() + "/analysis/all_functions/" + function + "_msgraph.txt");
	print_graphic(mgraph, out3); out3.close();
}

// main tester
/*
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "tot_info";
	//std::string function = "myfabs";
	TestType ttype = TestType::tot_info;

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

	// load MSG
	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto beg = cfiles.begin(), end = cfiles.end();
	while (beg != end) {
		// get code file
		const CodeFile & cfile = *(*(beg++));

		// test-module
		test_all_mutants_msg(root, cfile, funcs, cmutant, ctest, cscore);
		//test_sub_mutants_msg(root, cfile, function, funcs, cmutant, ctest, cscore);

		// end this file
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &ctrace;
	delete &cscore; delete &funcs;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	// end of all
	std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
*/