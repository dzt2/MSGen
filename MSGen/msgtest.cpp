#include "mclass.h"
#include "cfunc.h"
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
/* select operators */
static void select_operators(std::set<std::string> & operators) {
	operators.clear();

	/* E-selective set */
	operators.insert("u-OLLN");
	operators.insert("u-ORRN");
	operators.insert("u-OAAN");
	operators.insert("u-VDTR");
	operators.insert("u-VTWD");
	operators.insert("I-DirVarIncDec");
	operators.insert("I-IndVarIncDec");

	/* newly-subsuming ops */
	operators.insert("u-ORAN");
	operators.insert("u-OARN");
	operators.insert("u-CRCR");
	operators.insert("u-VLSR");
	operators.insert("u-VGSR");
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
/* print the information about MSG-structure */
static void print_ms_graph(const MS_Graph & graph, std::ostream & out) {
	// output structure summary
	int C = graph.size(), S = 0, E = 0;
	int M = graph.size_of_mutants();
	int minouts = 0, maxouts = 0;
	for (int k = 0; k < C; k++) {
		MSG_Node & node = graph.get_node(k);
		S += node.get_ou_port().degree();
		if (node.get_score_degree() == 0)
			E = node.get_mutants().number_of_mutants();
		else {
			minouts += (1 + node.get_ou_port().degree())
				* node.get_mutants().number_of_mutants();
			maxouts += (
				node.get_mutants().number_of_mutants() - 1
				+ node.get_ou_port().degree())
				* node.get_mutants().number_of_mutants();
		}
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
	int avg_min_outs = round(((double)minouts) / ((double)M));
	int avg_max_outs = round(((double)maxouts) / ((double)M));
	double min_outs_dens = ((double)avg_min_outs) / ((double)M);
	double max_outs_dens = ((double)avg_max_outs) / ((double)M);
	out << "\tMin-Outs: " << avg_min_outs << "/" << M << " (" << min_outs_dens << ")\n";
	out << "\tMax-Outs: " << avg_max_outs << "/" << M << " (" << max_outs_dens << ")\n";
	out << "-----------------------------\n\n";

	// output eof
	out << std::endl;
}

/* classify by operator */
void classify_by_operators(MutantSpace & space, std::map<Mutant::ID, MType> & typelib) {
	typelib.clear();
	Mutant::ID mid, n = space.number_of_mutants();
	for (mid = 0; mid < n; mid++) {
		Mutant & mutant = space.get_mutant(mid);
		typelib[mutant.get_id()] = mutant.get_operator();
	}
}
/* classify by functions */
void classify_by_functions(MutantSpace & space, const CFunctionSpace & funcs, std::map<Mutant::ID, MType> & typelib) {
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
/* classify by coverage */
void classify_by_coverage(MutantSpace & space, const MS_Graph & cgraph, std::map<Mutant::ID, MType> & typelib) {
	typelib.clear();
	Mutant::ID mid, n = space.number_of_mutants();
	for (mid = 0; mid < n; mid++) {
		if (cgraph.has_node_of(mid)) {
			MSG_Node & cnode = cgraph.get_node_of(mid);
			typelib[mid] = std::to_string(cnode.get_node_id());
		}
	}
}

/* trim the line-space */
static void trim_spaces(std::string & line) {
	std::string cache; char ch;
	for (int i = 0; i < line.length(); i++) {
		ch = line[i];
		if (ch != '\n') cache += ch;
	}
	line = cache;
}
/* whether x directly subsumes y */
static bool direct_subsume(MSG_Node & x, MSG_Node & y) {
	const MSG_Port & port = x.get_ou_port();
	for (int i = 0; i < port.degree(); i++) {
		MSG_Edge & edge = port.get_edge(i);
		MSG_Node & next = edge.get_target();
		if (&next == &y) return true;
	}
	return false;
}

/* classify graph */
static void classify_graph(const MS_Graph & source, MS_C_Graph & target, const std::map<Mutant::ID, MType> & typelib) {
	MS_C_Build builder;
	builder.open(typelib);
	builder.classify(source, target);
	builder.close();
}
/* print classify graph */
static void print_classify_graph(const MS_C_Graph & graph, std::ostream & out) {
	// for counters 
	size_t inner_eqs = 0, inter_eqs = 0;
	size_t inner_dss = 0, inter_dss = 0;
	size_t mutants = 0;

	MS_C_Node::ID cid, n = graph.size();
	for (cid = 0; cid < n; cid++) {
		MS_C_Node & x = graph.get_node(cid);
		size_t x_mutants = x.number_of_mutants();
		inner_eqs += x_mutants * (x_mutants - 1) ;
		mutants += x_mutants;

		const std::set<MS_C_Node *> & nexts = x.get_next_nodes();
		auto beg = nexts.begin(), end = nexts.end();
		while (beg != end) {
			MS_C_Node & y = *(*(beg++));
			size_t y_mutants = y.number_of_mutants();

			// equivalence between types 
			if (y.is_linked_to(x)) {
				if (x.get_type() == y.get_type()) {
					std::cout << "\t\t[error]Impossible case: " << x.get_type() << "\n";
					throw "error";
				}
				else inter_eqs += x_mutants * y_mutants;
			}
			// direct subsumption 
			else {
				if(x.get_type() == y.get_type())
					inner_dss ++;
				else inter_dss ++;
			}
		}
	}

	out << "\tC-Mutants: " << mutants << "\n";
	out << "\tC-Inner-EQ: " << inner_eqs << "/" << (inner_eqs + inter_eqs) << " (" << ((double)inner_eqs) / ((double)(inner_eqs + inter_eqs)) << ")\n";
	out << "\tC-Inter-EQ: " << inter_eqs << "/" << (inner_eqs + inter_eqs) << " (" << ((double)inter_eqs) / ((double)(inner_eqs + inter_eqs)) << ")\n";
	out << "\tC-Inner-DS: " << inner_dss << "/" << (inner_dss + inter_dss) << " (" << ((double)inner_dss) / ((double)(inner_dss + inter_dss)) << ")\n";
	out << "\tC-Inter-DS: " << inter_dss << "/" << (inner_dss + inter_dss) << " (" << ((double)inter_dss) / ((double)(inner_dss + inter_dss)) << ")\n";
	out << std::endl;
}
/* print inter-type equivalence */
static void print_mutants(const MS_Graph & graph, 
	const std::map<Mutant::ID, MType> & typelib, std::ostream & out) {
	out << "mid\tnode\toperator\ttype\tline\torigin\treplace\n";

	MutantSpace & space = graph.get_space();
	const TextBuild & txt = *(space.get_code_file().get_text());

	Mutant::ID n = space.number_of_mutants();
	for (Mutant::ID mid = 0; mid < n; mid++) {
		Mutant & mutant = space.get_mutant(mid);
		if (graph.has_node_of(mid)) {
			MSG_Node & node = graph.get_node_of(mid);
			if (node.get_score_degree() == 0) continue;
			auto iter = typelib.find(mid);

			const Mutation & mutation = mutant.
				get_mutation(mutant.get_orders() - 1);
			const CodeLocation & loc = mutation.get_location();
			std::string origin = loc.get_text_at();
			std::string replace = mutation.get_replacement();
			trim_spaces(origin); trim_spaces(replace);

			out << mid << "\t" << node.get_node_id() << "\t";
			out << mutant.get_operator() << "\t" << iter->second;
			out << "\t" << txt.lineOfIndex(loc.get_bias());
			out << "\t" << origin << "\t" << replace << "\n";
		}
	}

}
/* print mutant subsumption graph */
static void print_ms_edges(const MS_Graph & graph, std::ostream & out) {
	out << "source\ttarget\n";
	long mid, n = graph.size();
	for (mid = 0; mid < n; mid++) {
		MSG_Node & src = graph.get_node(mid);
		const MSG_Port & port = src.get_ou_port();
		for (int k = 0; k < port.degree(); k++) {
			MSG_Edge & edge = port.get_edge(k);
			MSG_Node & trg = edge.get_target();
			out << src.get_node_id() << 
				"\t" << trg.get_node_id() << "\n";
		}
	}
	out << std::endl;
}

/* test methods */
static void test_ms_graph(const File & root, const CodeFile & cfile, CFuncProject & funcs,
	CMutant & cmutant, CTest & ctest, CScore & cscore, MS_Graph & graph) {
	// get set of mutants and tests in project
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
	funcs.load_functions_for(cfile);
	CFunctionSpace & funcspace = funcs.get_function_space(cfile);

	// head information
	std::cout << "Load functions:\n\t{";
	const std::set<std::string> & funcnames = funcspace.get_function_names();
	auto beg = funcnames.begin(), end = funcnames.end();
	while (beg != end) {
		std::cout << *(beg++);
		if (beg != end) std::cout << "; ";
	}
	std::cout << "}\n";

	// get score vector producer | consumer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	FileScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	// MS-Graph-Build
	build_up_graph(graph, producer, consumer);

	// print graph
	// print_ms_graph(graph, std::cout);

	// output equivalence
	std::map<Mutant::ID, MType> typelib;
	classify_by_functions(mspace, funcspace, typelib);
	std::ofstream fout(root.get_path() + "/analysis/mutants.txt");
	print_mutants(graph, typelib, fout); fout.close();
	std::ofstream eout(root.get_path() + "/analysis/graph.txt");
	print_ms_edges(graph, eout); eout.close();

	// release resource
	ctest.delete_test_set(&tests); mspace.delete_set(&mutants);
}
static void test_ms_graph(const File & root, const CodeFile & cfile, CFuncProject & funcs, 
	CMutant & cmutant, CTest & ctest, CScore & cscore, const std::set<std::string> & operators, 
	MS_Graph & graph) {
	// get set of mutants and tests in project
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
	funcs.load_functions_for(cfile);
	CFunctionSpace & funcspace = funcs.get_function_space(cfile);

	// head information
	std::cout << "Load functions:\n\t{";
	const std::set<std::string> & funcnames = funcspace.get_function_names();
	auto beg = funcnames.begin(), end = funcnames.end();
	while (beg != end) {
		std::cout << *(beg++);
		if (beg != end) std::cout << "; ";
	}
	std::cout << "}\n";

	// get score vector producer | consumer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	FileScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	// select mutants
	std::set<Mutant::ID> mutset;
	select_mutants_by_operators(mspace, operators, mutset);
	ScoreFilter filter(producer, mutset);

	// MS-Graph-Build
	build_up_graph(graph, filter, consumer);

	// print graph
	// print_ms_graph(graph, std::cout);

	// output equivalence
	std::map<Mutant::ID, MType> typelib;
	classify_by_functions(mspace, funcspace, typelib);
	std::ofstream fout(root.get_path() + "/analysis/mutants.txt");
	print_mutants(graph, typelib, fout); fout.close();
	std::ofstream eout(root.get_path() + "/analysis/graph.txt");
	print_ms_edges(graph, eout); eout.close();

	// release resource
	ctest.delete_test_set(&tests); mspace.delete_set(&mutants);
}
static void test_cv_graph(const File & root, const CodeFile & cfile, CFuncProject & funcs,
	CMutant & cmutant, CTest & ctest, CTrace & ctrace, CScore & cscore, MS_Graph & graph) {
	// load coverage information
	CoverageSpace & covspac = ctrace.get_space();
	MutantSpace & mspace = cmutant.get_mutants_of(cfile);
	MutantSet & mutants = *(mspace.create_set()); mutants.complement();
	TestSet & tests = *(ctest.malloc_test_set()); tests.complement();
	funcs.load_functions_for(cfile);
	CFunctionSpace & funcspace = funcs.get_function_space(cfile);
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

	// outputs graph
	std::map<Mutant::ID, MType> typelib;
	classify_by_functions(mspace, funcspace, typelib);
	std::ofstream fout(root.get_path() + "/analysis/cmutants.txt");
	print_mutants(graph, typelib, fout); fout.close();
	std::ofstream eout(root.get_path() + "/analysis/cgraph.txt");
	print_ms_edges(graph, eout); eout.close();

	// release resource
	ctest.delete_test_set(&tests);
}
static void test_cv_ms_graphs(const MutantSpace & mspace, const MS_Graph & mgraph, const MS_Graph & cgraph) {
	Mutant::ID mid, n = mspace.number_of_mutants();
	std::cout << "/---------- validation ----------/\n";
	for (mid = 0; mid < n; mid++) {
		if (mgraph.has_node_of(mid)) {
			MSG_Node & snode = mgraph.get_node_of(mid);
			MSG_Node & cnode = cgraph.get_node_of(mid);
			if (!snode.get_score_vector().subsume(cnode.get_score_vector())) {
				std::cout << "\t(error): unable to understand mutant #" << mid << "\n";
			}
			//else std::cout << "\t(pass) mutant #" << mid << "\n";
		}
	}
	std::cout << "/------------------------------/\n";
}
static void test_cv_ms_equivalence(MutantSpace & mspace, const MS_Graph & mgraph, const MS_Graph & cgraph) {
	Mutant::ID i, j, n = mspace.number_of_mutants();
	size_t eqs = 0, subsuming = 0, subsumed = 0, nons = 0;
	size_t M = 0;
	for (i = 0; i < n; i++) {
		if (mgraph.has_node_of(i)) {
			M++;
			MSG_Node & ni = mgraph.get_node_of(i);
			MSG_Node & ci = cgraph.get_node_of(i);
			const BitSeq & bi = ci.get_score_vector();

			for (j = i + 1; j < n; j++) {
				if (mgraph.has_node_of(j)) {
					MSG_Node & nj = mgraph.get_node_of(j);
					MSG_Node & cj = cgraph.get_node_of(j);
					const BitSeq & bj = cj.get_score_vector();

					if (ni.get_node_id() == nj.get_node_id()) {
						if (bi.equals(bj)) eqs++;
						else if (bi.subsume(bj))
							subsuming++;
						else if (bj.subsume(bi))
							subsumed++;
						else nons++;
					}
				}
			}
		}
	}

	size_t total = eqs + subsuming + subsumed + nons;
	std::cout << "\tBlocks: " << cgraph.size() << "\n";
	std::cout << "\tMutants: " << M << "\n";
	std::cout << "\tEquivalent: " << eqs << "/" << total << " (" << ((double)eqs) / ((double)total) << ")\n";
	std::cout << "\tSubsumings: " << subsuming << "/" << total << " (" << ((double)subsuming) / ((double)total) << ")\n";
	std::cout << "\tSubsummeds: " << subsumed << "/" << total << " (" << ((double)subsumed) / ((double)total) << ")\n";
	std::cout << "\tNonRelated: " << nons << "/" << total << " (" << ((double)nons) / ((double)total) << ")\n";
}
static void evaluate_equivalence(MutantSpace & space, const MS_Graph & mgraph, const MS_Graph & cgraph) {
	std::map<Mutant::ID, MType> typelib;
	classify_by_coverage(space, cgraph, typelib);

	MS_C_Graph mc_graph;
	classify_graph(mgraph, mc_graph, typelib);

	size_t in_eqs = 0, po_eqs = 0, no_eqs = 0, dr_dss = 0;
	size_t in_dss = 0, po_dss = 0, ng_dss = 0, no_dss = 0;

	long n = mc_graph.size();
	for (long i = 0; i < n; i++) {
		/* get mutant mi and its coverage */
		MS_C_Node & ni = mc_graph.get_node(i);
		size_t mi = ni.number_of_mutants();
		unsigned ci = std::stoi(ni.get_type());
		MSG_Node & covi = cgraph.get_node(ci);
		const BitSeq & cbi = covi.get_score_vector();

		in_eqs += mi - 1;

		const std::set<MS_C_Node *> & nexts = ni.get_next_nodes();
		auto beg = nexts.begin(), end = nexts.end();
		while (beg != end) {
			/* get mutants mj and its coverage */
			MS_C_Node & nj = *(*(beg++));
			size_t mj = nj.number_of_mutants();
			unsigned cj = std::stoi(nj.get_type());
			MSG_Node & covj = cgraph.get_node(cj);
			const BitSeq & cbj = covj.get_score_vector();

			if (nj.is_linked_to(ni)) {	// equivalence
				if (ci == cj) {}
				else if (cbi.subsume(cbj)) 
					po_eqs++;
				else no_eqs++;
			}	
			else {						// direct subsumption
				if (ci == cj) in_dss++;
				else if (cbi.subsume(cbj)) {
					po_dss++;
					if (direct_subsume(covi, covj))
						dr_dss++;
				}
				else if (cbj.subsume(cbi))
					ng_dss++;
				else no_dss++;
			}
		}
	}

	std::cout << "\n---------- Relations ----------\n";
	std::cout << "\tIn_Eqs: " << in_eqs << "\n";
	std::cout << "\tPo_Eqs: " << po_eqs << "\n";
	std::cout << "\tNo_Eqs: " << no_eqs << "\n";
	std::cout << "\tIn_Dss: " << in_dss << "\n";
	std::cout << "\tDr_Dss: " << dr_dss << "\n";
	std::cout << "\tPo_Dss: " << po_dss << "\n";
	std::cout << "\tNg_Dss: " << ng_dss << "\n";
	std::cout << "\tNo_Dss: " << no_dss << "\n";
	std::cout << "-------------------------------\n";
}

// analysis methods
/* load the functions */
static void load_functions(CFuncProject & funcs, CMutant & cmutant, 
	const CodeFile & cfile, std::map<Mutant::ID, MType> & typelib) {
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
/* get the set of nodes in MSG that are subsumed by specific mutant */
static void get_subsumed_mutations(const MS_Graph & graph, Mutant::ID mid, std::set<MSG_Node *> & subsume_set) {
	subsume_set.clear();
	if (graph.has_node_of(mid)) {
		MSG_Node & node = graph.get_node_of(mid);
		if (node.get_score_degree() != 0) {
			subsume_set.insert(&node);
			std::queue<MSG_Node *> queue;

			queue.push(&node);
			while (!queue.empty()) {
				MSG_Node * next = queue.front(); queue.pop();
				const MSG_Port & port = next->get_ou_port();
				for (int k = 0; k < port.degree(); k++) {
					MSG_Edge & edge = port.get_edge(k);
					MSG_Node & trg = edge.get_target();
					if (subsume_set.count(&trg) == 0) {
						queue.push(&trg);
						subsume_set.insert(&trg);
					}
				}
			} // end while
		}
	}
}
/* evaluate the height of mutants among each block of CSG */
static void evaluate_height_blocks(const MS_Graph & mgraph, const MS_Graph & cgraph, std::map<MSG_Node *, double> & heights) {
	/* get the counts */
	std::map<MSG_Node *, unsigned> counters;
	MutantSpace & mspace = cgraph.get_space();
	Mutant::ID n = mspace.number_of_mutants();
	for (Mutant::ID mid = 0; mid < n; mid++) {
		if (cgraph.has_node_of(mid) && mgraph.has_node_of(mid)) {
			/* get the node of this mutant in coverage-graph */
			MSG_Node & ci = cgraph.get_node_of(mid);
			MSG_Node & mi = mgraph.get_node_of(mid);
			if (mi.get_score_degree() == 0) continue;
			if (ci.get_score_degree() == 0) continue;

			/* get the mutation height for this mutant */
			std::set<MSG_Node *> subsume_set;
			get_subsumed_mutations(mgraph, mid, subsume_set);

			/* count the height in block node */
			if (counters.count(&ci) == 0)
				counters[&ci] = 0;
			auto iter = counters.find(&ci);
			counters[&ci] = (iter->second) + subsume_set.size() - 1;
		}
	}

	/* calculate the average height */
	heights.clear();
	auto beg = counters.begin();
	auto end = counters.end();
	while (beg != end) {
		MSG_Node & ci = (*((beg)->first));
		unsigned count = (beg++)->second;
		double height = ((double)count) / ((double)
			ci.get_mutants().number_of_mutants());
		heights[&ci] = height;
	}
	counters.clear();
}
/* evaluate the number of mutants subsumed by those in the blocks */
static void evaluate_number_blocks(const MS_Graph & mgraph, const MS_Graph & cgraph, std::map<MSG_Node *, double> & numbers) {
	std::map<MSG_Node *, size_t> counters;
	MutantSpace & mspace = cgraph.get_space();
	Mutant::ID n = mspace.number_of_mutants();

	for (Mutant::ID mid = 0; mid < n; mid++) {
		if (cgraph.has_node_of(mid) && mgraph.has_node_of(mid)) {
			/* get the node of this mutant in coverage-graph */
			MSG_Node & ci = cgraph.get_node_of(mid);
			MSG_Node & mi = mgraph.get_node_of(mid);
			if (mi.get_score_degree() == 0) continue;
			if (ci.get_score_degree() == 0) continue;

			/* get the mutation height for this mutant */
			std::set<MSG_Node *> subsume_set;
			get_subsumed_mutations(mgraph, mid, subsume_set);

			/* get the number of mutants */
			size_t count = 0;
			auto beg = subsume_set.begin();
			auto end = subsume_set.end();
			while (beg != end) {
				MSG_Node & trg = *(*(beg++));
				count += trg.get_mutants().number_of_mutants();
			}

			/* count the height in block node */
			if (counters.count(&ci) == 0)
				counters[&ci] = 0;
			auto iter = counters.find(&ci);
			counters[&ci] = (iter->second) + count - 1;
		}
	}

	/* calculate average number of subsumed ones */
	numbers.clear();
	auto beg = counters.begin();
	auto end = counters.end();
	while (beg != end) {
		MSG_Node & ci = (*((beg)->first));
		unsigned count = (beg++)->second;
		double number = ((double)count) / ((double)
			ci.get_mutants().number_of_mutants());
		numbers[&ci] = number;
	}
	counters.clear();

}

/* std::cout analysis */
static void compare_height_coverage(std::map<MSG_Node *, double> & heights) {
	size_t pos = 0, neg = 0;

	auto beg = heights.begin();
	auto end = heights.end();
	while (beg != end) {
		MSG_Node & ni = *(beg->first); 
		double hi = (beg->second);
		const BitSeq & bi = ni.get_score_vector();

		auto beg2 = ++beg;
		while (beg2 != end) {
			MSG_Node & nj = *(beg2->first);
			double hj = (beg2->second);
			const BitSeq & bj = nj.get_score_vector();

			if (bi.subsume(bj)) {
				if (hi > hj) pos++;
				else neg++;
			}
			else if (bj.subsume(bi)) {
				if (hi > hj) neg++;
				else pos++;
			}
			beg2++;
		}
	}

	/* outputs */
	std::cout << "/---------- height analysis ----------/\n";
	std::cout << "\tPos-case: " << pos << "/" << (pos + neg) 
		<< "\t(" << ((double)pos) / ((double)(pos + neg)) << ")\n";
	std::cout << "\tNeg-case: " << neg << "/" << (pos + neg)
		<< "\t(" << ((double)neg) / ((double)(pos + neg)) << ")\n";
	std::cout << "/-------------------------------------/\n";
}

// outputs methods
/* print the information of graph */
static void print_graph(const MS_Graph & graph, std::ostream & out) {
	out << "node_id\tscore_degree\tmutants_size\tnext_set\n";

	size_t n = graph.size();
	for (size_t i = 0; i < n; i++) {
		MSG_Node & ni = graph.get_node(i);
		out << ni.get_node_id() << "\t";
		out << ni.get_score_degree() << "\t";
		out << ni.get_mutants().number_of_mutants() << "\t";

		const MSG_Port & port = ni.get_ou_port();
		for (int k = 0; k < port.degree(); k++) {
			MSG_Edge & edge = port.get_edge(k);
			MSG_Node & trg = edge.get_target();
			out << trg.get_node_id() << "; ";
		}
		out << "\n";
	}
}
/* output the height for each block in cgraph */
static void print_height_blocks(const std::map<MSG_Node *, double> & heights, std::ostream & out) {
	out << "block\tdegree\tmutants\theight\n";
	auto beg = heights.begin();
	auto end = heights.end();
	while (beg != end) {
		MSG_Node & block = *((beg)->first);
		double height = ((beg++)->second);
		out << block.get_node_id() << "\t";
		out << block.get_score_degree() << "\t";
		out << block.get_mutants().number_of_mutants() << "\t";
		out << height << "\n";
	}
}
/* output the number of mutant subsumed by a block  */
static void print_number_blocks(MutantSpace & mspace, const std::map<MSG_Node *, double> 
	& numbers, const std::map<MSG_Node *, double> & heights, std::ostream & out) {
	out << "block\tdegree\tmutants\tc-height\tnumbers\tpercentage\n";

	size_t total = mspace.number_of_mutants();
	double value = 0.0, values = 0.0, sum = 0.0;
	auto beg = numbers.begin(), end = numbers.end();
	while (beg != end) {
		MSG_Node * cnode = beg->first;
		size_t number = beg->second;
		value = ((double)number) / ((double)total);
		auto iter = heights.find(cnode);
		double height = iter->second;

		out << cnode->get_node_id() << "\t";
		out << cnode->get_score_degree() << "\t";
		out << cnode->get_mutants().number_of_mutants() << "\t";
		out << height << "\t";
		out << number << "\t";
		out << value << "\n";

		values += value * (cnode->get_mutants().number_of_mutants());
		sum += cnode->get_mutants().number_of_mutants(); beg++;
	}

	out << "\taverage-value: " << values << "/" << sum << " (" << values / sum << ")\n";
}

/* main method */
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

	// load MSG
	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto beg = cfiles.begin(), end = cfiles.end();
	while (beg != end) {
		// get the next file of code
		const CodeFile & cfile = *(*(beg++));
		std::map<Mutant::ID, MType> funclib;
		load_functions(funcs, cmutant, cfile, funclib);
		std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

		// analysis declarations
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		MS_Graph mgraph(mspace), cgraph(mspace);

		// mutant subsumption graph construction 
		load_ms_graph(root, cfile, funcs, cmutant, ctest, cscore, mgraph);
		load_cs_graph(root, cfile, cmutant, ctest, ctrace, cscore, cgraph);
		std::ofstream msm_out(root.get_path() + "/analysis/msg_mutants.txt");
		std::ofstream msg_out(root.get_path() + "/analysis/msg_graph.txt");
		print_mutants(mgraph, funclib, msm_out); msm_out.close();
		print_graph(mgraph, msg_out); msg_out.close();
		std::ofstream csm_out(root.get_path() + "/analysis/csg_mutants.txt");
		std::ofstream csg_out(root.get_path() + "/analysis/csg_graph.txt");
		print_mutants(cgraph, funclib, csm_out); csm_out.close();
		print_graph(cgraph, csg_out); csg_out.close();

		// height analysis based on MSG and CSG
		std::map<MSG_Node *, double> heights;
		evaluate_height_blocks(mgraph, cgraph, heights);
		std::ofstream chigh(root.get_path() + "/analysis/cov_height.txt");
		print_height_blocks(heights, chigh); chigh.close(); 
		std::map<MSG_Node *, double> numbers;
		std::map<MSG_Node *, double> cheight;
		evaluate_height_blocks(cgraph, cgraph, cheight);
		evaluate_number_blocks(mgraph, cgraph, numbers);
		std::ofstream cnumb(root.get_path() + "/analysis/cov_number.txt");
		print_number_blocks(mspace, numbers, cheight, cnumb); cnumb.close();

		// std-out analysis 
		compare_height_coverage(heights);

		// correlations
		// test subsumption graph for mutants
		// test_ms_graph(root, cfile, funcs, cmutant, ctest, cscore, mgraph);
		// test subsumption graph for coverage
		// test_cv_graph(root, cfile, funcs, cmutant, ctest, ctrace, cscore, cgraph);
		// validate
		// test_cv_ms_graphs(mspace, mgraph, cgraph);
		// print-relations
		// evaluate_equivalence(mspace, mgraph, cgraph);

		// end this file
		heights.clear(); cheight.clear(); numbers.clear();
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &ctrace;
	delete &cscore; delete &funcs;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
