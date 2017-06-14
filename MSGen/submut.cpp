#include "submut.h"

TypedMutantSet::TypedMutantSet(MSGraph & g) : graph(g),
	mutants(g.get_class_set().get_mutants()) {
	
	MutantSpace & mspace = mutants.get_space();
	stubborn_mutants = mspace.create_set();
	hard_mutants = mspace.create_set();
	subsuming_mutants = mspace.create_set();
	subsumed_mutants = mspace.create_set();

	stubborn_cluster = nullptr;
	hard_clusters.clear();
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

	const MuHierarchy & hierarchy = graph.get_hierarchy();
	size_t degree = hierarchy.get_degrees()[0];
	if (stubborn_cluster != nullptr) degree = hierarchy.get_degrees()[1];
	const std::set<MuCluster *> & H1 = hierarchy.get_clusters_of(degree);
	auto hbeg = H1.begin(), hend = H1.end();
	while (hbeg != hend) {
		MuCluster * cluster = *(hbeg++);
		hard_clusters.insert(cluster);
	}

	Mutant::ID mid = 0, mnum = mspace.number_of_mutants();
	while (mid < mnum) {
		if (graph.has_cluster_of(mid)) {
			MuCluster & cluster = graph.get_cluster_of(mid);

			/* stubborn mutant */
			if (&cluster == stubborn_cluster) {
				stubborn_mutants->add_mutant(mid);
			}
			/* hard mutants in H[1] (kulled by least tests) */
			else if (hard_clusters.count(&cluster) > 0) {
				hard_mutants->add_mutant(mid);
				subsuming_mutants->add_mutant(mid);
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
	mspace.delete_set(hard_mutants);

	subsuming_clusters.clear();
	hard_clusters.clear();
	subsumed_clusters.clear();
}

void TypedOutputter::open(const File & root) {
	close();	/* close the outputter */

	/* iterate the root directory */
	if (root.is_directory()) {
		const std::vector<File *> & files = root.list_files();
		for (int i = 0; i < files.size(); i++) {
			File & file = *(files[i]);
			if (file.is_directory()) {
				const std::string & name
					= file.get_local_name();
				if (name == "analysis") {
					dir = &file; break;
				}
			}
		}
	}

	/* validation */
	if (dir == nullptr) {
		CError error(CErrorType::InvalidArguments, 
			"TypedOutputter::open", "no directory ../analysis is found");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
}
void TypedOutputter::output_mutants(const TypedMutantSet & tmutants) {
	std::ofstream out1(dir->get_path() + "/stubborn_mutants.txt");
	output_mutants(tmutants.get_stubborn_mutants(), 
		tmutants.get_graph(), out1); out1.close();

	std::ofstream out2(dir->get_path() + "/hard_mutants.txt");
	output_mutants(tmutants.get_hard_mutants(), 
		tmutants.get_graph(), out2); out2.close();

	std::ofstream out3(dir->get_path() + "/subsuming_mutants.txt");
	output_mutants(tmutants.get_subsuming_mutants(), 
		tmutants.get_graph(), out3); out3.close();
}
void TypedOutputter::output_distribution(const TypedMutantSet & tmutants) {
	std::ofstream out1(dir->get_path() + "/mutant_operator.txt");
	output_distribute_operator(tmutants, out1); out1.close();

	std::ofstream out2(dir->get_path() + "/mutant_location.txt");
	output_distribute_location(tmutants, out2); out2.close();

	std::ofstream out3(dir->get_path() + "/mutant_operator_location.txt");
	output_distribute_operator_location(tmutants, out3); out3.close();
}
void TypedOutputter::output_mutants(const MutantSet & mutants, 
	/* declarations */
	const MSGraph & graph, std::ostream & out) {
	MutantSpace & mspace = mutants.get_space();
	const TextBuild & text = *(mspace.get_code_file().get_text());

	/* iterate all the mutants in the set */
	Mutant::ID mid, num = mspace.number_of_mutants();
	for (mid = 0; mid < num; mid++) {
		if (mutants.has_mutant(mid)) {
			/* get the next mutants */
			Mutant & mutant = mspace.get_mutant(mid);
			const std::string & oprt = mutant.get_operator();
			MuCluster & cluster = graph.get_cluster_of(mid);

			/* mutant-summary */
			out << "Mutant (" << mutant.get_id() << ")\n";
			out << "\toperator: " << oprt << "\n";
			out << "\tcluster: " << cluster.get_id() << "\n";
			out << "\tdegree: " << cluster.get_score_degree() << "\n";

			/* mutations for mutant */
			for (int i = 0; i < mutant.get_orders(); i++) {
				const Mutation & mutation = mutant.get_mutation(i);
				const CodeLocation & location = mutation.get_location();

				out << "\tmutation[" << i << "]:\n";
				out << "\t\t-location: " << 
					text.lineOfIndex(location.get_bias()) << " line\n";
				out << "\t\t-original: \"" << location.get_text_at() << "\"\n";
				out << "\t\t-replace: \"" << mutation.get_replacement() << "\"\n";
			}

			out << "\n";	/* end line */
		}
	}
}
void TypedOutputter::output_distribute_operator(
	const TypedMutantSet & tmutants, std::ostream & out) {
	/* declarations */
	const MutantSet & all_mutants = tmutants.get_mutants();
	const MutantSet & stubborn_mutants = tmutants.get_stubborn_mutants();
	const MutantSet & subsuming_mutants = tmutants.get_subsuming_mutants();
	std::map<std::string, int *> oprt_numbers;

	const MutantSpace & mspace = all_mutants.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	for (mid = 0; mid < num; mid++) {
		/* get next mutant */
		Mutant & mutant = mspace.get_mutant(mid);
		const std::string & oprt = mutant.get_operator();

		/* get the numbers for counting */
		int * numbers;
		if (oprt_numbers.count(oprt) == 0) {
			numbers = new int[3];
			numbers[0] = 0; 
			numbers[1] = 0;
			numbers[2] = 0; 
			oprt_numbers[oprt] = numbers;
		}
		else {
			auto iter = oprt_numbers.find(oprt);
			numbers = iter->second;
		}

		/* count the numbers */
		if (all_mutants.has_mutant(mid)) {
			/* stubborn mutants */
			if (stubborn_mutants.has_mutant(mid))
				numbers[0]++;
			/* subsuming mutants */
			else if (subsuming_mutants.has_mutant(mid))
				numbers[1]++;
			/* subsumed mutants */
			else numbers[2]++;	
		}
	}

	/* output operators */
	out << "operator\tStubborn\tSubsuming\tSubsumed\tTotal\n";
	auto beg = oprt_numbers.begin(), end = oprt_numbers.end();
	while (beg != end) {
		const std::string & oprt = beg->first;
		int * numbers = beg->second; beg++;

		out << oprt << "\t" << numbers[0] << "\t";
		out << numbers[1] << "\t" << numbers[2] << "\t";
		out << numbers[0] + numbers[1] + numbers[2] << "\n";

		delete numbers;
	}
	out << std::endl;
}
void TypedOutputter::output_distribute_location(
	const TypedMutantSet & tmutants, std::ostream & out) {
	/* declarations */
	const MutantSet & all_mutants = tmutants.get_mutants();
	const MutantSet & stubborn_mutants = tmutants.get_stubborn_mutants();
	const MutantSet & subsuming_mutants = tmutants.get_subsuming_mutants();
	std::map<std::string, int *> key_numbers;
	std::map<std::string, const CodeLocation *> key_locations;

	/* count the numbers */
	const MutantSpace & mspace = all_mutants.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	for (mid = 0; mid < num; mid++) {
		/* get next mutant and its location-key */
		Mutant & mutant = mspace.get_mutant(mid);
		const Mutation & mutation = 
			mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & location = mutation.get_location();
		std::string key = std::to_string(location.get_bias());
		key += ":" + std::to_string(location.get_length());

		/* get numbers for counting */
		int * numbers;
		if (key_numbers.count(key) == 0) {
			numbers = new int[3];
			numbers[0] = 0;
			numbers[1] = 0;
			numbers[2] = 0;
			key_numbers[key] = numbers;
			key_locations[key] = &location;
		}
		else {
			auto iter = key_numbers.find(key);
			numbers = iter->second;
		}

		/* couting for mutants */
		if (all_mutants.has_mutant(mid)) {
			/* for stubborn mutants */
			if (stubborn_mutants.has_mutant(mid))
				numbers[0]++;
			/* for subsuming mutants */
			else if (subsuming_mutants.has_mutant(mid))
				numbers[1]++;
			/* for subsumed mutants */
			else numbers[2]++;
		}
	} /* end for */

	/* output lines */
	const TextBuild & text = *(mspace.get_code_file().get_text());
	out << "bias\tlength\tline\tstubborn\tsubsuming\tsubsumed\ttotal\n";
	auto beg = key_numbers.begin(), end = key_numbers.end();
	while (beg != end) {
		/* get next line */
		const std::string & key = beg->first;
		int * numbers = beg->second; beg++;
		auto iter = key_locations.find(key);
		const CodeLocation & loc = *(iter->second);

		/* output line */
		out << loc.get_bias() << "\t" << loc.get_length() << "\t";
		out << text.lineOfIndex(loc.get_bias()) << "\t" << numbers[0];
		out << "\t" << numbers[1] << "\t" << numbers[2] << "\t";
		out << numbers[0] + numbers[1] + numbers[2] << "\n";

		delete numbers;
	}
	out << std::endl;
}
void TypedOutputter::output_distribute_operator_location(
	const TypedMutantSet & tmutants, std::ostream & out) {
	/* declarations */
	const MutantSet & all_mutants = tmutants.get_mutants();
	const MutantSet & stubborn_mutants = tmutants.get_stubborn_mutants();
	const MutantSet & subsuming_mutants = tmutants.get_subsuming_mutants();
	std::map<std::string, int *> key_numbers;
	std::map<std::string, std::string> key_operators;
	std::map<std::string, const CodeLocation *> key_locations;

	/* count the numbers */
	const MutantSpace & mspace = all_mutants.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	for (mid = 0; mid < num; mid++) {
		/* get next mutant and operator-location-key */
		Mutant & mutant = mspace.get_mutant(mid);
		const Mutation & mutation =
			mutant.get_mutation(mutant.get_orders() - 1);
		const std::string oprt = mutant.get_operator();
		const CodeLocation & location = mutation.get_location();
		
		/* get the key for this mutant */
		std::string key = std::to_string(location.get_bias());
		key += ":" + std::to_string(location.get_length());
		key += ":" + oprt;

		/* get numbers for counting */
		int * numbers;
		if (key_numbers.count(key) == 0) {
			numbers = new int[3];
			numbers[0] = 0;
			numbers[1] = 0;
			numbers[2] = 0;
			key_numbers[key] = numbers;
			key_operators[key] = oprt;
			key_locations[key] = &location;
		}
		else {
			auto iter = key_numbers.find(key);
			numbers = iter->second;
		}

		/* couting for mutants */
		if (all_mutants.has_mutant(mid)) {
			/* for stubborn mutants */
			if (stubborn_mutants.has_mutant(mid))
				numbers[0]++;
			/* for subsuming mutants */
			else if (subsuming_mutants.has_mutant(mid))
				numbers[1]++;
			/* for subsumed mutants */
			else numbers[2]++;
		}
	} /* end for */

	/* output */
	const TextBuild & text = *(mspace.get_code_file().get_text());
	out << "operator\tbias\tlength\tline\tstubborn\tsubsuming\tsubsumed\ttotal\n";
	auto beg = key_numbers.begin(), end = key_numbers.end();
	while (beg != end) {
		/* get next line */
		const std::string & key = beg->first;
		int * numbers = beg->second; beg++;
		auto iter = key_locations.find(key);
		const CodeLocation & loc = *(iter->second);
		const std::string & oprt =
			(key_operators.find(key)->second);

		/* output line */
		out << oprt << "\t" << loc.get_bias() << "\t" << loc.get_length() << "\t";
		out << text.lineOfIndex(loc.get_bias()) << "\t" << numbers[0] << "\t";
		out << numbers[1] << "\t" << numbers[2] << "\t";
		out << numbers[0] + numbers[1] + numbers[2] << "\n";

		delete numbers;	/* release numbers */
	}
	out << std::endl;
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
	out << " \tStubborn\tSubsuming\tHard\tSubsumed\tTotal\n";
	out << "mutants\t" << tmutants.get_stubborn_mutants().number_of_mutants();
	out << "\t" << tmutants.get_subsuming_mutants().number_of_mutants();
	out << "\t" << tmutants.get_hard_mutants().number_of_mutants();
	out << "\t" << tmutants.get_subsumed_mutants().number_of_mutants();
	out << "\t" << tmutants.get_mutants().number_of_mutants() << "\n";

	out << "clusters\t";
	tmutants.get_stubborn_cluster();
	out << "1\t";
	out << tmutants.get_subsuming_clusters().size() << "\t";
	out << tmutants.get_hard_clusters().size() << "\t";
	out << tmutants.get_subsumed_clusters().size() << "\t";
	out << tmutants.get_graph().size() << "\n";
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

/* test main method */
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/"; 
	std::string prname = "minmax"; TestType ttype = TestType::general;

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
		TypedMutantSet tmutants(graph);

		// output MSG
		std::ofstream out(prefix + prname + "/statistics.txt");
		out << "Mutant Subsumption Graph\n";
		printMSG(graph, out);
		efficiencyMSG(graph, out);
		printTMS(tmutants, out);
		out.close();

		// create TypedMutantSet
		TypedOutputter tout; tout.open(root);
		tout.output_mutants(tmutants);
		tout.output_distribution(tmutants);
		tout.close();
	}

	// delete memory
	delete &cmutant; delete &ctest; 
	delete &program; delete &root;

	// exit
	std::cout << "Press any key to exit...\n"; getchar(); exit(0);
}