#include "sgraph.h"
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

// analysis methods
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

/* calculator methods */
/* get the set of nodes in MSG|CSG that are subsumed by target node (not including itself) */
static bool calculate_nodes_subsumed_by(MSG_Node & root, std::set<MSG_Node *> & ans) {
	ans.clear();
	if (root.get_score_degree() == 0) 
		return false;
	else {
		std::queue<MSG_Node *> qlist;
		qlist.push(&root);
		while (!qlist.empty()) {
			MSG_Node & src = *(qlist.front()); qlist.pop();
			const MSG_Port & port = src.get_ou_port();
			for (int k = 0; k < port.degree(); k++) {
				MSG_Edge & edge = port.get_edge(k);
				MSG_Node & trg = edge.get_target();

				if (ans.count(&trg) == 0) {
					ans.insert(&trg);
					qlist.push(&trg);
				}
			}
		}
		return true;
	}
}
/* calculate the set of nodes in MSG|CSG that subsume the leaf node */
static bool calculate_nodes_subsummings(MSG_Node & leaf, std::set<MSG_Node *> & ans) {
	ans.clear();
	if (leaf.get_score_degree() == 0)
		return false;
	else {
		std::queue<MSG_Node *> qlist;
		qlist.push(&leaf);
		while (!qlist.empty()) {
			MSG_Node & trg = *(qlist.front()); qlist.pop();
			const MSG_Port & port = trg.get_in_port();
			for (int k = 0; k < port.degree(); k++) {
				MSG_Edge & edge = port.get_edge(k);
				MSG_Node & src = edge.get_source();

				if (src.get_score_degree() > 0
					&& ans.count(&src) == 0) {
					ans.insert(&src);
					qlist.push(&src);
				}
			}
		}
		return true;
	}
}
/* calculate the number of mutants in the nodes set */
static size_t calculate_mutants_of_nodes(const std::set<MSG_Node *> & nodes) {
	size_t sum = 0;

	auto beg = nodes.begin();
	auto end = nodes.end();
	while (beg != end) {
		MSG_Node & node = *(*(beg++));
		sum += node.get_mutants().number_of_mutants();
	}

	return sum;
}
/* calculate the number of nodes and mutants subsumed (or subsuming) the target one */
static void calculate_nodes_and_mutants(MSG_Node & node, size_t ans[4]) {
	std::set<MSG_Node *> subsummeds;
	std::set<MSG_Node *> subsumings;
	calculate_nodes_subsummings(node, subsumings);
	calculate_nodes_subsumed_by(node, subsummeds);
	ans[0] = subsummeds.size();
	ans[1] = calculate_mutants_of_nodes(subsummeds);
	ans[2] = subsumings.size();
	ans[3] = calculate_mutants_of_nodes(subsumings);
}
/* calculate the set of subsuming mutants */
static void calculate_subsuming_nodes(MS_Graph & graph, std::set<MSG_Node *> & subsumings) {
	// get equivalent ones
	subsumings.clear();
	size_t k, n = graph.size(), e;
	for (e = n, k = 0; k < n; k++) {
		MSG_Node & node = graph.get_node(k);
		const MSG_Port & port = node.get_in_port();

		if (node.get_score_degree() == 0) {
			e = k; break;
		}
		else if (port.degree() == 0) {
			subsumings.insert(&node);
		}
	}

	// compute subsumings 
	if (e != n) {
		MSG_Node & root = graph.get_node(e);
		const MSG_Port & port = root.get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & next = edge.get_target();
			subsumings.insert(&next);
		}
	}
}
/* get the average utility of each mutant in set */
static double calculate_mutants_utility(MS_Graph & graph, const std::set<Mutant::ID> & mutants) {
	std::set<MSG_Node *> nodes;
	auto beg = mutants.begin();
	auto end = mutants.end();

	double utility = 0.0;
	while (beg != end) {
		Mutant::ID mid = *(beg++);
		if (graph.has_node_of(mid)) {
			MSG_Node & node = graph.get_node_of(mid);
			std::set<MSG_Node *> subsumings;
			std::set<MSG_Node *> subsummeds;
			calculate_nodes_subsumed_by(node, subsummeds);
			calculate_nodes_subsummings(node, subsumings);
			size_t A = subsummeds.size();
			size_t B = subsumings.size();
			if (A == 0 && B == 0) continue;
			double use = ((double)A) / ((double)(A + B));
			utility = utility + use;
		}
	}

	return utility / mutants.size();
}

/* file-output methods */
/* print the information of mutants in the graph */
static void print_mutants(const MS_Graph & msg, const MS_Graph & csg, 
	const std::set<MSG_Node *> & subsumings, std::ostream & out) {
	MutantSpace & mspace = msg.get_space();
	Mutant::ID n = mspace.number_of_mutants();
	const CodeFile & cfile = mspace.get_code_file();
	const TextBuild & text = *(cfile.get_text());

	out << "id\ttype\toperator\tline\torigin\treplace\tmsg_node\tcsg_node\n";
	for (Mutant::ID mid = 0; mid < n; mid++) {
		if (msg.has_node_of(mid) && csg.has_node_of(mid)) {
			MSG_Node & mnode = msg.get_node_of(mid);
			MSG_Node & cnode = csg.get_node_of(mid);
			if (mnode.get_score_degree() != 0) {
				Mutant & mutant = mspace.get_mutant(mid);
				const std::string & op = mutant.get_operator();
				const Mutation & mutation = mutant.
					get_mutation(mutant.get_orders() - 1);
				const CodeLocation & loc = mutation.get_location();
				std::string origin = loc.get_text_at();
				std::string replace = mutation.get_replacement();
				trim_spaces(origin); trim_spaces(replace);

				out << mid << "\t";

				if (mnode.get_score_degree() == 0)
					out << "equivalent\t";
				else if (subsumings.count(&mnode) > 0)
					out << "subsuming\t";
				else out << "subsummed\t";

				out << op << "\t";
				out << text.lineOfIndex(loc.get_bias()) << "\t";
				out << origin << "\t" << replace << "\t";
				out << mnode.get_node_id() << "\t";
				out << cnode.get_node_id() << "\n";
			}
		}
	}
	out << std::endl;
}
/* print the information of nodes in graph */
static void print_graph(const MS_Graph & graph, std::ostream & out) {
	std::set<MSG_Node *> subsumings;
	std::set<MSG_Node *> subsummeds;
	long n = graph.size();

	out << "node\tmutants\tdegree\tsubsuming-nodes\tsubsuming-mutants\tsubsummed-nodes\tsubsummed-mutants\tnext_set\n";
	for (long id = 0; id < n; id++) {
		MSG_Node & node = graph.get_node(id);
		if (node.get_score_degree() > 0) {
			out << id << "\t";
			out << node.get_mutants().number_of_mutants() << "\t";
			out << node.get_score_degree() << "\t";
			
			calculate_nodes_subsumed_by(node, subsummeds);
			calculate_nodes_subsummings(node, subsumings);
			size_t subsummeds_num = calculate_mutants_of_nodes(subsummeds);
			size_t subsumings_num = calculate_mutants_of_nodes(subsumings);
			out << subsumings.size() << "\t";
			out << subsumings_num << "\t";
			out << subsummeds.size() << "\t";
			out << subsummeds_num << "\t";

			const MSG_Port & port = node.get_ou_port();
			for (int k = 0; k < port.degree(); k++) {
				MSG_Edge & edge = port.get_edge(k);
				MSG_Node & trg = edge.get_target();
				out << trg.get_node_id() << "; ";
			}

			out << "\n";
		}
	}
	out << std::endl;
}
/* print {degree-utility-subsumed} */
static void print_block(MS_Graph & csg, MS_Graph & msg, std::ostream & out) {
	/* title */
	out << "blockid\tcover-height\tcover-mutants\tcover-utility\t";
	out << "msg-avg-degree\tmsg-avg-subsummeds\tmsg-avg-utility\t";
	out << "mut-avg-degree\tmut-avg-subsummeds\tmut-avg-utility\n";

	/* get the arguments for each MSG-node */
	std::map<MSG_Node *, size_t *> msg_args;
	size_t msg_n = msg.size();
	for (size_t i = 0; i < msg_n; i++) {
		MSG_Node & mnode = msg.get_node(i);
		size_t * arg = new size_t[4];
		calculate_nodes_and_mutants(mnode, arg);
		msg_args[&mnode] = arg;
	}

	/* relate csg to msg */
	MSG_Relation relations(csg, msg);

	/* get each pair and calculate their arguments */
	size_t csg_n = csg.size(); size_t args[4];
	for (size_t i = 0; i < csg_n; i++) {
		MSG_Node & cnode = csg.get_node(i);
		if (relations.has_related_targets(cnode)) {
			/* blockid */
			out << cnode.get_node_id() << "\t";
			
			/* coverage-argument */
			calculate_nodes_and_mutants(cnode, args);
			out << args[0] << "\t" << args[1] << '\t';
			out << ((double)args[0]) / ((double)(args[0] + args[2])) << "\t";

			/* get pairs */
			const std::set<MSG_Pair *> & pairs = relations.get_related_targets(cnode);
			double CN = pairs.size(), MN = cnode.get_mutants().number_of_mutants();

			/* calculate mutant arguments */
			double AD = 0.0, AS = 0.0, AU = 0.0;
			double MD = 0.0, MS = 0.0, MU = 0.0;
			double CD = 0.0, CS = 0.0, CU = 0.0;
			auto beg = pairs.begin(), end = pairs.end();
			while (beg != end) {
				MSG_Pair & pair = *(*(beg++));
				double PN = pair.get_mutants().number_of_mutants();
				MSG_Node & mnode = pair.get_target();
				if (msg_args.count(&mnode) > 0) {
					size_t * margs = (msg_args.find(&mnode))->second;

					/* current arguments */
					CD = mnode.get_score_degree(); CS = margs[1];
					CU = ((double)margs[0]) / ((double)(margs[0] + margs[2]));

					/* calculate the average-arguments */
					AD += CD, AS += CS, AU += CU;
					MD += CD * PN / MN; 
					MS += CS * PN / MN;
					MU += CU * PN / MN;
				}
			} // end while

			/* calculate average arguments */
			AD = AD / CN; AS = AS / CN; AU = AU / CN;

			/* output arguments */
			out << AD << "\t" << AS << "\t" << AU << "\t";
			out << MD << "\t" << MS << "\t" << MU << "\n";
		} // end if
	} // end for 
	out << std::endl;

	/* delete the arguments for MSG nodes */
	auto mbeg = msg_args.begin();
	auto mend = msg_args.end();
	while (mbeg != mend)
		delete (*(mbeg++)).second;
	msg_args.clear();
}
/* print {value; nodes; mutants;}*/
static void print_value(MS_Graph & msg, std::ostream & out) {
	/* declarations */
	std::map<MSG_Node *, size_t> node_value;
	size_t n = msg.size(), k = 0, value;
	std::set<MSG_Node *> subsummeds;

	/* compute the value for each node */
	for (k = 0; k < n; k++) {
		/* get the next node from MSG */
		MSG_Node & node = msg.get_node(k);
		if (node.get_score_degree() == 0) continue;

		/* calculate the value */
		calculate_nodes_subsumed_by(node, subsummeds);
		value = calculate_mutants_of_nodes(subsummeds);
		value += node.get_mutants().number_of_mutants() - 1;

		/* record the value */
		node_value[&node] = value;
	}
	
	/* count values */
	std::map<size_t, size_t> value_nodes;
	std::map<size_t, size_t> value_mutants;
	auto beg = node_value.begin();
	auto end = node_value.end();
	while (beg != end) {
		/* get next node and value */
		MSG_Node & node = *(beg->first);
		size_t value = beg->second;

		/* insert */
		if (value_nodes.count(value) == 0) 
			value_nodes[value] = 0;
		if (value_mutants.count(value) == 0)
			value_mutants[value] = 0;

		auto iter1 = value_nodes.find(value);
		value_nodes[value] = (iter1->second) + 1;
		auto iter2 = value_mutants.find(value);
		value_mutants[value] = (iter2->second) 
			+ node.get_mutants().number_of_mutants();

		beg++;
	}
	node_value.clear();

	/* output information */
	out << "value\tnodes\tmutants\n";
	auto beg1 = value_nodes.begin();
	auto end1 = value_nodes.end();
	while (beg1 != end1) {
		// get the value of 
		size_t value = beg1->first;
		size_t nodes = beg1->second;
		auto iter3 = value_mutants.find(value);
		size_t mutants = iter3->second;

		out << value << "\t" << nodes << 
			"\t" << mutants << "\n";

		beg1++;
	}
	out << std::endl;

	/* clear */ 
	value_nodes.clear();
	value_mutants.clear();
}

/* analysis methods */
static void evaluate_csg_nodes(MS_Graph & csg, MS_Graph & msg, std::map<MSG_Node *, double *> & ans) {
	/* evaluate node in MSG */
	std::map<MSG_Node *, size_t *> msglib;
	size_t msgn = msg.size(), csgn = csg.size();
	for (size_t i = 0; i < msgn; i++) {
		MSG_Node & node = msg.get_node(i);
		size_t * arg = new size_t[4];
		calculate_nodes_and_mutants(node, arg);
		msglib[&node] = arg;
	}

	/* build up relations from CSG to MSG */
	MSG_Relation relations(csg, msg);

	/* calculate the arguments for node in CSG */
	for (size_t i = 0; i < csgn; i++) {
		MSG_Node & cnode = csg.get_node(i);
		if (relations.has_related_targets(cnode)) {
			/* get pairs */
			const std::set<MSG_Pair *> & pairs = relations.get_related_targets(cnode);
			double CN = pairs.size(), MN = cnode.get_mutants().number_of_mutants();

			/* calculate mutant arguments */
			double AD = 0.0, AS = 0.0, AU = 0.0;
			double MD = 0.0, MS = 0.0, MU = 0.0;
			double CD = 0.0, CS = 0.0, CU = 0.0;
			auto beg = pairs.begin(), end = pairs.end();
			while (beg != end) {
				MSG_Pair & pair = *(*(beg++));
				double PN = pair.get_mutants().number_of_mutants();
				MSG_Node & mnode = pair.get_target();
				if (msglib.count(&mnode) > 0) {
					size_t * margs = (msglib.find(&mnode))->second;

					/* current arguments */
					CD = mnode.get_score_degree(); CS = margs[1];
					CU = ((double)margs[0]) / ((double)(margs[0] + margs[2]));

					/* calculate the average-arguments */
					AD += CD, AS += CS, AU += CU;
					MD += CD * PN / MN;
					MS += CS * PN / MN;
					MU += CU * PN / MN;
				}
			} // end while

			  /* calculate average arguments */
			AD = AD / CN; AS = AS / CN; AU = AU / CN;

			/* record the results */
			double * res = new double[4];
			res[0] = AS, res[1] = AU;
			res[2] = MS, res[3] = MU;
			ans[&cnode] = res;
		}
	} // end for

	  /* delete the arguments for MSG nodes */
	auto mbeg = msglib.begin();
	auto mend = msglib.end();
	while (mbeg != mend)
		delete (*(mbeg++)).second;
	msglib.clear();
}
/* output the positive-negative proportions */
static void ouput_block(MS_Graph & csg, MS_Graph & msg) {
	std::map<MSG_Node *, double *> csglib;
	evaluate_csg_nodes(csg, msg, csglib);

	size_t AS_Pos = 0, AS_Neg = 0;
	size_t AU_Pos = 0, AU_Neg = 0;
	size_t MS_Pos = 0, MS_Neg = 0;
	size_t MU_Pos = 0, MU_Neg = 0;

	auto beg = csglib.begin();
	auto end = csglib.end();
	while (beg != end) {
		MSG_Node & src = *((beg++)->first);
		const BitSeq & sbits = src.get_score_vector();
		if (csglib.count(&src) == 0) continue;
		auto iter1 = csglib.find(&src);
		double * src_args = iter1->second;

		auto beg2 = beg;
		while (beg2 != end) {
			MSG_Node & trg = *((beg2++)->first);
			const BitSeq & tbits = trg.get_score_vector();
			if (csglib.count(&trg) == 0) continue;
			auto iter2 = csglib.find(&trg);
			double * trg_args = iter2->second;

			if (sbits.subsume(tbits)) {
				/* AS-Positive-Negative */
				if (src_args[0] >= trg_args[0]) 
					AS_Pos++;
				else AS_Neg++;
				/* AU-Positive-Negative */
				if (src_args[1] >= trg_args[1])
					AU_Pos++;
				else AU_Neg++;
				/* MS-Positive-Negative */
				if (src_args[2] >= trg_args[2])
					MS_Pos++;
				else MS_Neg++;
				/* MU-Positive-Negative */
				if (src_args[3] >= trg_args[3])
					MU_Pos++;
				else MU_Neg++;
			}
			else if (tbits.subsume(sbits)) {
				/* AS-Positive-Negative */
				if (src_args[0] > trg_args[0])
					AS_Neg++;
				else AS_Pos++;
				/* AU-Positive-Negative */
				if (src_args[1] > trg_args[1])
					AU_Neg++;
				else AU_Pos++;
				/* MS-Positive-Negative */
				if (src_args[2] > trg_args[2])
					MS_Neg++;
				else MS_Pos++;
				/* MU-Positive-Negative */
				if (src_args[3] > trg_args[3])
					MU_Neg++;
				else MU_Pos++;
			}
		}
	}

	/* delete resources */
	beg = csglib.begin();
	end = csglib.end();
	while (beg != end)
		delete ((beg++)->second);
	csglib.clear();

	/* output */
	std::cout << "/-------------------------------------------/\n";
	std::cout << "\tMSG-AVG-Sub: " << AS_Pos << "\t" << AS_Neg << "\n";
	std::cout << "\tMSG-AVG-Uty: " << AU_Pos << "\t" << AU_Neg << "\n";
	std::cout << "\tMut-AVG-Sub: " << MS_Pos << "\t" << MS_Neg << "\n";
	std::cout << "\tMut-AVG-Uty: " << MU_Pos << "\t" << MU_Neg << "\n";
	std::cout << "/------------------------------------------/\n";
}
/* indistinguish(inner|dominated|others) direct(inner|dominated|others) */
static void output_subsumption(MS_Graph & msg, MS_Graph & csg) {
	/* counters */
	size_t eq_inner = 0, eq_domat = 0, eq_other = 0;
	size_t di_inner = 0, di_domat = 0, di_other = 0;

	/* construct relations */
	MSG_Relation relations(msg, csg);
	const std::map<std::string, MSG_Pair *> & pairs = relations.get_pairs();

	/* indistinguishable relations */
	size_t msg_n = msg.size();
	for (size_t i = 0; i < msg_n; i++) {
		MSG_Node & node = msg.get_node(i);
		if (relations.has_related_targets(node)) {
			const std::set<MSG_Pair *> & pairs 
				= relations.get_related_targets(node);
			auto beg1 = pairs.begin(), end = pairs.end();
			while (beg1 != end) {
				MSG_Pair & pair1 = *(*(beg1++));
				size_t mutants = pair1.get_mutants().number_of_mutants();
				//eq_inner = eq_inner + mutants * (mutants - 1) / 2;
				eq_inner = eq_inner + mutants - 1;

				MSG_Node & B1 = pair1.get_target();
				const BitSeq & bit1 = B1.get_score_vector();

				for (auto beg2 = beg1; beg2 != end; beg2++) {
					MSG_Pair & pair2 = *(*(beg2));
					MSG_Node & B2 = pair2.get_target();
					const BitSeq & bit2 = B2.get_score_vector();

					if (bit1.subsume(bit2)) eq_domat++;
					else if (bit2.subsume(bit1)) eq_domat++;
					else eq_other++;
				}

				
			}
		}
	}
	
	/* direct subsumption */
	for (size_t i = 0; i < msg_n; i++) {
		MSG_Node & node = msg.get_node(i);
		if (node.get_score_degree() == 0) continue;
		else if (relations.has_related_targets(node)) {
			const std::set<MSG_Pair *> & src_blocks
				= relations.get_related_targets(node);

			const MSG_Port & port = node.get_ou_port();
			for (int k = 0; k < port.degree(); k++) {
				MSG_Edge & edge = port.get_edge(k);
				MSG_Node & next = edge.get_target();
				if (!relations.has_related_targets(next)) continue;

				const std::set<MSG_Pair *> & trg_blocks
					= relations.get_related_targets(next);

				auto src_beg = src_blocks.begin();
				auto src_end = src_blocks.end();
				while (src_beg != src_end) {
					MSG_Pair & src_pair = *(*(src_beg++));
					size_t src_mutants = src_pair.
						get_mutants().number_of_mutants();
					MSG_Node & B1 = src_pair.get_target();
					const BitSeq & C1 = B1.get_score_vector();

					auto trg_beg = trg_blocks.begin();
					auto trg_end = trg_blocks.end();
					while (trg_beg != trg_end) {
						MSG_Pair & trg_pair = *(*(trg_beg++));
						size_t trg_mutants = trg_pair.
							get_mutants().number_of_mutants();
						MSG_Node & B2 = trg_pair.get_target();
						const BitSeq & C2 = B2.get_score_vector();

						if (&B1 == &B2)
							di_inner++;
						else if (C1.subsume(C2))
							di_domat++;
						else
							di_other++;
					}
				}
			}
		}
	}

	/* output conclusions */
	std::cout << "\n/---------- subsumption distribute ----------/\n";
	std::cout << "\tEquivalence-inner: " << eq_inner << "\n";
	std::cout << "\tEquivalence-domat: " << eq_domat << "\n";
	std::cout << "\tEquivalence-other: " << eq_other << "\n";
	std::cout << "\tDirectSubsm-inner: " << di_inner << "\n";
	std::cout << "\tDirectSubsm-domat: " << di_domat << "\n";
	std::cout << "\tDirectSubsm-other: " << di_other << "\n";
	std::cout << "/---------------------------------------------/\n";

}

static void select_mutants(std::set<Mutant::ID> & mutants, File & root) {
	std::ifstream in(root.get_path() + "/analysis/selects.txt");
	std::string line;
	while (std::getline(in, line)) {
		if (line.empty()) continue;
		else {
			Mutant::ID mid = std::atoi(line.c_str());
			mutants.insert(mid);
		}
	}
	in.close();
}
static void select_mutants(size_t s, size_t n, std::set<Mutant::ID> & mutants) {
	mutants.clear();
	while (s-- != 0) {
		Mutant::ID mid;
		do {
			/* get standard seed */
			double sed = std::rand();
			sed = sed / RAND_MAX;

			/* get mutant id */
			mid = std::floor(sed * n);
		} while (mutants.count(mid) > 0);
		mutants.insert(mid);
	}
}
static void output_score(MS_Graph & mgraph, CTest & ctest, 
	const std::set<Mutant::ID> & mutants) {
	TestSet * tests = ctest.malloc_test_set();

	// compute minimal test set
	MSG_Tester tester;
	tester.open(mgraph);
	tester.gen_tests(mutants, *tests);

	// compute tests mutation score
	double score = tester.eval_score(*tests);
	double domsc = tester.eval_dom_score(*tests);

	// utility
	double utility = calculate_mutants_utility(mgraph, mutants);

	// precision & recall
	std::set<MSG_Node *> roots, recalled;
	calculate_subsuming_nodes(mgraph, roots);
	double precision = 0, recall = 0;
	auto beg = mutants.begin(), end = mutants.end();
	while (beg != end) {
		Mutant::ID mid = *(beg++);
		if (mgraph.has_node_of(mid)) {
			MSG_Node & node = mgraph.get_node_of(mid);
			long node_id = node.get_node_id();
			if (roots.count(&node) > 0) {
				recalled.insert(&node);
				precision = precision + 1.0;
			}
		}
	}
	precision = precision / mutants.size();
	recall = ((double)recalled.size()) / ((double)roots.size());

	// output
	std::cout << "\n/--------- Score Analysis ----------/\n";
	std::cout << "\tMutants:\t" << mutants.size() << "\n";
	std::cout << "\tTests:  \t" << tests->size() << "\n";
	std::cout << "\tScore:  \t" << score << "\n";
	std::cout << "\tDScore: \t" << domsc << "\n";
	std::cout << "\tPrecis: \t" << precision << "\n";
	std::cout << "\tRecall: \t" << recall << "\n";
	std::cout << "\tUtility:\t" << utility << "\n";
	std::cout << "/-------------------------------------/\n";

	// delete
	ctest.delete_test_set(tests); tester.close();
}

/* main method */
/*
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "bubble";
	TestType ttype = TestType::general;

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
		std::map<Mutant::ID, std::string> funclib;
		load_functions(funcs, cmutant, cfile, funclib);
		std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

		// analysis declarations
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		MS_Graph mgraph(mspace), cgraph(mspace);

		// mutant subsumption graph construction 
		load_ms_graph(root, cfile, funcs, cmutant, ctest, cscore, mgraph);
		load_cs_graph(root, cfile, cmutant, ctest, ctrace, cscore, cgraph);

		// calculate
		std::set<MSG_Node *> subsumings;
		calculate_subsuming_nodes(mgraph, subsumings);

		// output information about mutants
		std::ofstream out1(root.get_path() + "/analysis/mutants.txt");
		print_mutants(mgraph, cgraph, subsumings, out1); out1.close();
		std::ofstream out2(root.get_path() + "/analysis/msgraph.txt");
		print_graph(mgraph, out2); out2.close();
		std::ofstream out3(root.get_path() + "/analysis/csgraph.txt");
		print_graph(cgraph, out3); out3.close();
		std::ofstream out4(root.get_path() + "/analysis/mblocks.txt");
		print_block(cgraph, mgraph, out4); out4.close();
		std::ofstream out5(root.get_path() + "/analysis/mvalues.txt");
		print_value(mgraph, out5); out5.close();

		// stdout analysis
		ouput_block(cgraph, mgraph);
		output_subsumption(mgraph, cgraph);

		// user-select mutants
		std::set<Mutant::ID> mutants;
		select_mutants(mutants, root);
		output_score(mgraph, ctest, mutants);
		// rand-select mutants
		size_t n = mspace.number_of_mutants();
		//size_t s = mutants.size();
		size_t s = n * 0.1;
		select_mutants(s, n, mutants);
		output_score(mgraph, ctest, mutants);

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
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &ctrace;
	delete &cscore; delete &funcs;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
*/