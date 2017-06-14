#include "submut.h"

TypedMutantSet::TypedMutantSet(MSGraph & g) : graph(g),
	mutants(g.get_class_set().get_mutants()) {
	
	MutantSpace & mspace = mutants.get_space();
	stubborn_mutants = mspace.create_set();
	subsuming_mutants = mspace.create_set();
	subsumed_mutants = mspace.create_set();

	stubborn_cluster = nullptr;
	subsuming_clusters.clear();
	subsumed_clusters.clear();

	MuCluster::ID cid = 0, num = graph.size();
	while (cid < num) {
		MuCluster & cluster = graph.get_cluster(cid++);

		if (cluster.get_score_degree() == 0) {
			if (stubborn_cluster != nullptr) {
				CError error(CErrorType::Runtime, 
					"TypedMutantSet::TypedMutantSet", 
					"Duplicated cluster: all-zeros");
				CErrorConsumer::consume(error);
				exit(CErrorType::Runtime);
			}
			else stubborn_cluster = &cluster;
		}
		else subsumed_clusters.insert(&cluster);
	}

	if (stubborn_cluster == nullptr) {
		const std::set<MuCluster *> & roots = graph.get_roots();
		auto beg = roots.begin(), end = roots.end();
		while (beg != end) {
			subsuming_clusters.insert(*(beg));
			subsumed_clusters.erase(*(beg));
			beg++;
		}
	}
	else {
		const std::vector<MuSubsume> & edges 
			= stubborn_cluster->get_ou_port().get_edges();
		auto beg = edges.begin(), end = edges.end();
		while (beg != end) {
			const MuSubsume & edge = *(beg++);
			MuCluster & trg = edge.get_target();
			subsuming_clusters.insert(&trg);
			subsumed_clusters.erase(&trg);
		}
	}

	Mutant::ID mid = 0, mnum = mspace.number_of_mutants();
	while (mid < mnum) {
		if (graph.has_cluster_of(mid)) {
			MuCluster & cluster = graph.get_cluster_of(mid);

			/* stubborn mutant */
			if (&cluster == stubborn_cluster) {
				stubborn_mutants->add_mutant(mid);
			}
			/* subsuming mutants */
			else if (subsuming_clusters.count(&cluster) > 0) {
				subsuming_mutants->add_mutant(mid);
			}
			/* subsumed mutants */
			else if (subsumed_clusters.count(&cluster) > 0) {
				subsumed_mutants->add_mutant(mid);
			}
			/* invalid classification */
			else {
				CError error(CErrorType::Runtime,
					"TypedMutantSet::TypedMutantSet",
					"Internal error occurs.");
				CErrorConsumer::consume(error);
				exit(CErrorType::Runtime);
			}
		}
		mid++;	/* roll to next mutant */
	}

}
TypedMutantSet::~TypedMutantSet() {
	MutantSpace & mspace = graph.get_class_set().get_mutants().get_space();
	mspace.delete_set(stubborn_mutants);
	mspace.delete_set(subsuming_mutants);
	mspace.delete_set(subsumed_mutants);
}

/* APIs for project models */
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
/* construct MSG for the specified mutants */
static void constructMSG(MSGraph & graph, MutantSet & mutants, TestSet & tests) {
	// get models 
	const CMutant & cmutant = mutants.get_space().get_project();
	const CTest & ctest = tests.get_space().get_project();
	const File & root = cmutant.get_code_space().get_project().get_root();
	const CodeFile & cfile = mutants.get_space().get_code_file();

	// get score function
	CScore & cscore = *(new CScore(root, cmutant, ctest));
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	ScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	// create unlinker MSG
	MSGBuilder builder;
	builder.install(graph);
	builder.build_up(producer, consumer);
	builder.uninstall();

	// link the nodes in MSG
	MSGLinker linker;
	linker.connect(graph, MSGLinker::top_down);
}

/* output APIs */
static void printMSG(const MSGraph & graph, std::ostream & out) {
	size_t edges = 0;
	MuCluster::ID cid = 0, num = graph.size();
	while (cid < num) {
		MuCluster & cluster = graph.get_cluster(cid++);
		edges += cluster.get_ou_port().get_degree();
	}

	out << "Total Summary of Graph\n";
	out << "\t(1) Clusters: \t" << graph.size() << "\n";
	out << "\t(2) Hierarchy:\t" << graph.get_hierarchy().size_of_degress() << "\n";
	out << "\t(3) Subsumes: \t" << edges << "\n";
}
static void printTMS(TypedMutantSet & tmutants, std::ostream & out) {
	out << " \tStubborn\tSubsuming\tSubsumed\tTotal\n";
	out << "mutants\t" << tmutants.get_stubborn_mutants().number_of_mutants();
	out << "\t" << tmutants.get_subsuming_mutants().number_of_mutants();
	out << "\t" << tmutants.get_subsumed_mutants().number_of_mutants();
	out << "\t" << tmutants.get_mutants().number_of_mutants() << "\n";

	out << "clusters\t";
	tmutants.get_stubborn_cluster();
	out << "1\t";
	out << tmutants.get_subsuming_clusters().size() << "\t";
	out << tmutants.get_subsumed_clusters().size() << "\t";
	out << tmutants.get_graph().size() << "\n";
}
static void mutantsByOperator(const MutantSet & mutants, std::ostream & out) {
	MuClassSet class_set(mutants);

	MuClassifierByOperator classifier;
	classifier.classify(class_set);
	
	out << "\tMutant Classification by Operators\n";
	const std::map<MuFeature, MuClass *> & 
		classes = class_set.get_classes();
	auto beg = classes.begin();
	auto end = classes.end();
	while (beg != end) {
		MuClass & _class = *((beg++)->second);
		std::string * oprt = (std::string *) _class.get_feature();

		out << "\t" << *oprt << ": \t" << 
			_class.get_mutants().number_of_mutants() << "\n";
	}
	out << std::endl;
}
static void mutantsByLocation(const MutantSet & mutants, std::ostream & out) {
	MuClassSet class_set(mutants);

	MuClassifierByLocation classifier;
	classifier.classify(class_set);

	out << "\tMutant Classification by Operators\n";
	const std::map<MuFeature, MuClass *> &
		classes = class_set.get_classes();
	auto beg = classes.begin();
	auto end = classes.end();
	while (beg != end) {
		MuClass & _class = *((beg++)->second);
		CodeLocation * loc = (CodeLocation *) _class.get_feature();
		out << "\t[" << loc->get_bias() << ", " << loc->get_length() << "]: \t" <<
			_class.get_mutants().number_of_mutants() << "\n";
	}
	out << std::endl;
}
static void efficiencyMSG(const MSGraph & graph, std::ostream & out) {
	out << "Efficiency Analysis for Mutant Subsumption Graph\n";

	size_t C = graph.size();
	size_t M = graph.get_class_set().get_mutants().number_of_mutants();
	size_t Sm = M * (M - 1) / 2, Sc = C * (C - 1) / 2, Sh, Sa, Se;
	
	const MuHierarchy & hierarchy = graph.get_hierarchy();
	int hlength = hierarchy.size_of_degress() - 1; 
	size_t cnt = 0, level_num; Sh = 0;
	while (hlength >= 0) {
		level_num = hierarchy.get_clusters_at(hlength--).size();
		Sh += level_num * cnt; cnt += level_num;
	}

	Sa = times;

	MuCluster::ID cid = 0, cnum = graph.size(); Se = 0;
	while (cid < cnum) {
		MuCluster & ci = graph.get_cluster(cid++);
		Se += ci.get_ou_port().get_degree();
	}

	out << "\tS[mutants] = \t" << Sm << "\n";
	out << "\tS[cluster] = \t" << Sc << "\n";
	out << "\tS[hierarchy] = \t" << Sh << "\n";
	out << "\tS[algorithm] = \t" << Sa << "\n";
	out << "\tS[minimal] = \t" << Se << "\n";
}

/* analysis APIs */
static void stubborn_operator(const TypedMutantSet & tmutants, std::ostream & out) {
	std::map<std::string, int *> oprt_map;
	const MutantSpace & mspace = tmutants.get_mutants().get_space();
	Mutant::ID mid = 0, num = mspace.number_of_mutants();

	const MutantSet & all_mutants = tmutants.get_mutants();
	const MutantSet & stubborn_mutants = tmutants.get_stubborn_mutants();
	const MutantSet & subsuming_mutants = tmutants.get_subsuming_mutants();

	while (mid < num) {
		Mutant & mutant = mspace.get_mutant(mid++);
		const std::string & oprt = mutant.get_operator();

		int * numbers;
		if (oprt_map.count(oprt) == 0) {
			numbers = new int[3];
			numbers[0] = 0; 
			numbers[1] = 0;
			numbers[2] = 0;
			oprt_map[oprt] = numbers;
		}
		else {
			auto iter = oprt_map.find(oprt);
			numbers = iter->second;
		}

		if (all_mutants.has_mutant(mid - 1)) {
			if (stubborn_mutants.has_mutant(mid - 1)) 
				numbers[0]++;
			else if (subsuming_mutants.has_mutant(mid - 1))
				numbers[1]++;
			else numbers[2]++;
		}
	}

	auto beg = oprt_map.begin(), end = oprt_map.end();
	out << "Operators\tStubborn\tSubsuming\tSubsumed\tTotals\n";
	while (beg != end) {
		const std::string & oprt = (beg)->first;
		int * numbers = (beg++)->second;
		out << oprt << "\t" << numbers[0] << "\t" << numbers[1] << "\t" 
			<< numbers[2] << "\t" << (numbers[0] + numbers[1] + numbers[2]) << "\n";
		delete numbers;
	}
	oprt_map.clear();
}
static void stubborn_location(const TypedMutantSet & tmutants, std::ostream & out) {
	std::map<std::string, int *> loc_mutants;
	std::map<std::string, const CodeLocation *> name_loc;
	MutantSpace & mspace = tmutants.get_mutants().get_space();
	Mutant::ID mid = 0, num = mspace.number_of_mutants();

	const MutantSet & all_mutants = tmutants.get_mutants();
	const MutantSet & stubborn_mutants = tmutants.get_stubborn_mutants();
	const MutantSet & subsuming_mutants = tmutants.get_subsuming_mutants();

	/* count the numbers and locations for mutants */
	for (mid = 0; mid < num; mid++) {
		/* get next mutant */
		Mutant & mutant = mspace.get_mutant(mid);

		/* extract its location and name */
		const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & location = mutation.get_location();
		std::string name = std::to_string(location.get_bias());
		name += ":" + std::to_string(location.get_length());

		/* record location and derive its numbers */
		int * numbers;
		if (name_loc.count(name) == 0) {
			numbers = new int[3];
			numbers[0] = 0;
			numbers[1] = 0;
			numbers[2] = 0;
			loc_mutants[name] = numbers;
			name_loc[name] = &location;
		}
		else {
			auto iter = loc_mutants.find(name);
			numbers = iter->second;
		}

		/* count the numbers */
		if (all_mutants.has_mutant(mid)) {
			if (stubborn_mutants.has_mutant(mid))
				numbers[0]++;
			else if (subsuming_mutants.has_mutant(mid))
				numbers[1]++;
			else numbers[2]++;
		}
	}

	/* get text */
	const TextBuild & text = *(mspace.get_code_file().get_text());

	/* output */
	out << "Bias\tLength\tLine\tStubborn\tSubsuming\tSubsumed\tTotal\n";
	auto beg = name_loc.begin(), end = name_loc.end();
	while (beg != end) {
		/* get next line */
		const std::string & name = beg->first;
		const CodeLocation & loc = *((beg++)->second);
		auto iter = loc_mutants.find(name);
		int * numbers = iter->second;

		/* print line */
		out << loc.get_bias() << "\t" << loc.get_length() << "\t";
		out << text.lineOfIndex(loc.get_bias()) << "\t";
		out << numbers[0] << "\t" << numbers[1] << "\t";
		out << numbers[2] << "\t" << numbers[0] + numbers[1] + numbers[2] << "\n";

		/* release memory */ delete numbers;
	}
	out << std::endl;
}
static void stubborn_location_operator(const TypedMutantSet & tmutants, std::ostream & out) {
	MutantSpace & mspace = tmutants.get_mutants().get_space();
	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	const TextBuild & text = *(mspace.get_code_file().get_text());

	const MutantSet & all_mutants = tmutants.get_mutants();
	const MutantSet & stubborn_mutants = tmutants.get_stubborn_mutants();
	const MutantSet & subsuming_mutants = tmutants.get_subsuming_mutants();

	std::map<std::string, std::string> name_oprt;
	std::map<std::string, const CodeLocation *> name_loc;
	std::map<std::string, int *> name_numbers;

	for (mid = 0; mid < num; mid++) {
		/* get next mutant, its location, operator and key */
		Mutant & mutant = mspace.get_mutant(mid);
		const std::string & oprt = mutant.get_operator();
		const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & location = mutation.get_location();
		std::string key = oprt + ":" + std::to_string(location.get_bias());
		key += ":" + std::to_string(location.get_length());

		/* get numbers */
		int * numbers;
		if (name_numbers.count(key) == 0) {
			numbers = new int[3];
			numbers[0] = 0;
			numbers[1] = 0;
			numbers[2] = 0;
			name_numbers[key] = numbers;
			name_oprt[key] = oprt;
			name_loc[key] = &location;
		}
		else {
			auto iter = name_numbers.find(key);
			numbers = iter->second;
		}

		/* count the numbers */
		if (all_mutants.has_mutant(mid)) {
			if (stubborn_mutants.has_mutant(mid))
				numbers[0]++;
			else if (subsuming_mutants.has_mutant(mid))
				numbers[1]++;
			else numbers[2]++;
		}
	}

	/* output */
	out << "Operator\tBias\tLength\tLine\tStubborn\tSubsuming\tSubsumed\tTotal\n";
	auto beg = name_numbers.begin(), end = name_numbers.end();
	while (beg != end) {
		const std::string & key = beg->first;
		int * numbers = beg->second; beg++;
		const std::string & oprt 
			= (name_oprt.find(key))->second;
		const CodeLocation & location
			= *(name_loc.find(key)->second);

		out << oprt << "\t" << location.get_bias() << "\t";
		out << location.get_length() << "\t";
		out << text.lineOfIndex(location.get_bias()) << "\t";
		out << numbers[0] << "\t" << numbers[1] << "\t";
		out << numbers[2] << "\t" << numbers[0] + numbers[1] + numbers[2] << "\n";

		delete numbers;
	}

	out << std::endl;
}
static void stubborn_mutants(const TypedMutantSet & tmutants, std::ostream & out) {
	const MutantSet & stubborn_mutants = 
		tmutants.get_stubborn_mutants();
	MutantSpace & mspace = stubborn_mutants.get_space();
	MSGraph & graph = tmutants.get_graph();

	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	const TextBuild & text = *(mspace.get_code_file().get_text());

	for (mid = 0; mid < num; mid++) {
		if (stubborn_mutants.has_mutant(mid)) {
			Mutant & mutant = mspace.get_mutant(mid);
			const std::string & oprt = mutant.get_operator();
			
			out << "Mutant {" << mutant.get_id() << "}\n";
			out << "\tOperator: " << oprt << "\n";
			out << "\tCluster: " << graph.get_cluster_of(mid).get_id() << "\n";
			for (int i = 0; i < mutant.get_orders(); i++) {
				const Mutation & mutation = mutant.get_mutation(i);
				const CodeLocation & location = mutation.get_location();
				const std::string & origin_code = location.get_text_at();
				const std::string & replace_code = mutation.get_replacement();

				out << "\tMutation[" << i << "]:\n";
				out << "\t\tLocation: " << text.lineOfIndex(location.get_bias()) << " line\n";
				out << "\t\tOriginal: \"" << origin_code << "\"\n";
				out << "\t\tReplace: \"" << replace_code << "\"\n";
			}
			out << "\n";
		}
	}
	out << std::endl;
}
static void subsuming_mutants(const TypedMutantSet & tmutants, std::ostream & out) {
	const MutantSet & subsuming_mutants =
		tmutants.get_subsuming_mutants();
	MutantSpace & mspace = subsuming_mutants.get_space();
	MSGraph & graph = tmutants.get_graph();

	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	const TextBuild & text = *(mspace.get_code_file().get_text());

	for (mid = 0; mid < num; mid++) {
		if (subsuming_mutants.has_mutant(mid)) {
			Mutant & mutant = mspace.get_mutant(mid);
			const std::string & oprt = mutant.get_operator();

			out << "Mutant {" << mutant.get_id() << "}\n";
			out << "\tOperator: " << oprt << "\n";
			out << "\tCluster: " << graph.get_cluster_of(mid).get_id() << "\n";
			for (int i = 0; i < mutant.get_orders(); i++) {
				const Mutation & mutation = mutant.get_mutation(i);
				const CodeLocation & location = mutation.get_location();
				const std::string & origin_code = location.get_text_at();
				const std::string & replace_code = mutation.get_replacement();

				out << "\tMutation[" << i << "]:\n";
				out << "\t\tLocation: " << text.lineOfIndex(location.get_bias()) << " line\n";
				out << "\t\tOriginal: \"" << origin_code << "\"\n";
				out << "\t\tReplace: \"" << replace_code << "\"\n";
			}
			out << "\n";
		}
	}
	out << std::endl;
}

/* test main method */
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/"; 
	std::string prname = "tcas"; 
	TestType ttype = TestType::tcas;

	// create code-project, mutant-project, test-project
	File & root = *(new File(prefix + prname));
	CProgram & program = *(new CProgram(root));
	CTest & ctest = *(new CTest(ttype, root, program.get_exec()));
	CMutant & cmutant = *(new CMutant(root, program.get_source()));

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
		MutantSet & mutants = *(mspace.create_set());
		TestSet & tests = *(ctest.malloc_test_set());
		mutants.complement(); tests.complement();

		std::cout << "Load file: \"" << cfile.get_file().get_path() << "\"\n";

		// create MSG
		MSGraph graph(mutants);
		constructMSG(graph, mutants, tests);

		// output MSG
		std::ofstream out(prefix + prname + "/statistics.txt");
		out << "Mutant Subsumption Graph\n";
		printMSG(graph, out);
		efficiencyMSG(graph, out);

		// create TypedMutantSet
		TypedMutantSet tmutants(graph);
		out << "\nTyped Mutants Set:\n";
		printTMS(tmutants, out);
		out << std::endl;

		// print classification
		out << "Classification for subsuming mutants\n";
		mutantsByOperator(tmutants.get_subsuming_mutants(), out);
		mutantsByLocation(tmutants.get_subsuming_mutants(), out);

		/* end output */ out.close();

		// try to output cluster-operator
		std::ofstream oprt_out(prefix + prname + "/mutant_operator.txt");
		stubborn_operator(tmutants, oprt_out); oprt_out.close();
		std::ofstream loc_out(prefix + prname + "/mutant_location.txt");
		stubborn_location(tmutants, loc_out); loc_out.close();
		std::ofstream oprt_loc_out(prefix + prname + "/mutant_operator_location.txt");
		stubborn_location_operator(tmutants, oprt_loc_out); oprt_loc_out.close();

		// output stubborn/subsuming mutants 
		std::ofstream stub_out(prefix + prname + "/stubborn.txt");
		stubborn_mutants(tmutants, stub_out); stub_out.close();
		std::ofstream subs_out(prefix + prname + "/subsuming.txt");
		subsuming_mutants(tmutants, subs_out); subs_out.close();
	}

	// delete memory
	delete &cmutant; delete &ctest; 
	delete &program; delete &root;

	// exit
	std::cout << "Press any key to exit...\n"; getchar(); exit(0);
}