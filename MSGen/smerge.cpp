#include "smerge.h"

void SuMutantSet::update_mutants() {
	/* initialize the mutant set */
	mutants->clear(); clusters.clear();

	/* get the nodes in graph that subsume all the others (may be equivalent) */
	const std::set<MuCluster *> & roots = graph.get_roots();
	auto beg = roots.begin(), end = roots.end();

	/* compute the equivalent cluster (if) */
	MuCluster * eqs = nullptr;
	if (roots.size() == 1) {
		MuCluster * cluster = *(beg++);
		if (cluster->get_score_degree() == 0)
			eqs = cluster;
	}
	equivalent_cluster = eqs;
	
	/* nodes directly subsumed by equivalent cluster are subsuming */
	if (eqs != nullptr) {
		const std::vector<MuSubsume> & edges 
			= eqs->get_ou_port().get_edges();
		auto ebeg = edges.begin(), eend = edges.end();
		while (ebeg != eend) {
			const MuSubsume & edge = *(ebeg++);
			MuCluster & target = edge.get_target();
			const MutantSet & tmuts = target.get_mutants();
			mutants->disjunct(tmuts); clusters.insert(&target);
		}
	}
	/* otherwise, the roots in graph are subsuming mutants */
	else {
		beg = roots.begin(), end = roots.end();
		while (beg != end) {
			MuCluster * cluster = *(beg++);
			const MutantSet & tmuts = cluster->get_mutants();
			mutants->disjunct(tmuts); clusters.insert(cluster);
		}
	}
}

void SuMutantMerger::open(SuMutantSet & result) {
	answer = &result;

	MSGraph & graph 
		= result.get_graph();
	builder.open(graph);
}
void SuMutantMerger::append(SuMutantSet & smuts) {
	/* get the subsuming mutants and its graph */
	const MutantSet & mutants = smuts.get_subsuming_mutants();
	Mutant::ID mid = 0, num = mutants.number_of_mutants();
	
	MSGraph & graph = smuts.get_graph();
	while (mid < num) {
		if (mutants.has_mutant(mid)) {
			MuCluster & cluster = graph.get_cluster_of(mid);
			const BitSeq & bits = cluster.get_score_vector();
			builder.add(mid, bits);
		}
		/* to the next mutant */ mid = mid + 1;
	}
}
void SuMutantMerger::extract() {
	builder.link(); answer->update_mutants();
}

const SuMutantSet & SuMutantExperimentCore::get_subsuming_mutants(const std::string & key) const {
	if (mut_map.count(key) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"SuMutantExperimentCore::get_subsuming_mutants", 
			"Invalid key: \"" + key + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = mut_map.find(key);
		SuMutantSet * smut = iter->second;
		return *smut;
	}
}
SuMutantSet & SuMutantExperimentCore::get_mutants(const std::string & key) {
	if (mut_map.count(key) == 0) {
		MSGraph * graph = new MSGraph(mspace);
		SuMutantSet * smut = new SuMutantSet(*graph);
		mut_map[key] = smut; keys.insert(key);
		return *smut;
	}
	else {
		auto iter = mut_map.find(key);
		SuMutantSet * smut = iter->second;
		return *smut;
	}
}
void SuMutantExperimentCore::clear_mutants() {
	auto beg = mut_map.begin();
	auto end = mut_map.end();
	while (beg != end) {
		SuMutantSet * smut = (beg++)->second;
		MSGraph & graph = smut->get_graph();
		delete smut; delete &graph;
	}
	mut_map.clear(); keys.clear();
}

void SuMutantExperimentDriver::derive_operator_I(
	ScoreProducer & producer, ScoreConsumer & consumer) {
	/* builders for initializing MSG according to operator */
	std::map<std::string, MSGBuilder *> builders;
	MutantSpace & mspace = core->get_space();

	/* create ulinker graph for mutants of each operator */
	ScoreVector * score_vector;
	while ((score_vector = producer.produce()) != nullptr) {
		/* get the next mutant */
		const BitSeq & bits = score_vector->get_vector();
		Mutant::ID mid = score_vector->get_mutant();
		Mutant & mutant = mspace.get_mutant(mid);
		const std::string & oprt = mutant.get_operator();

		/* get the subsuming set for the operator */
		SuMutantSet & smut = core->get_mutants(oprt);
		MSGraph & graph = smut.get_graph();

		/* get the builder */
		MSGBuilder * builder;
		if (builders.count(oprt) == 0) {
			builder = new MSGBuilder();
			builder->open(graph);
			builders[oprt] = builder;
		}
		else {
			auto iter = builders.find(oprt);
			builder = iter->second;
		}

		/* append the mutant into the graph-builder */
		builder->add(mid, bits);

		/* consume the vector */
		consumer.consume(score_vector);
	} /* end while */

	/* construct linked graph */
	auto beg = builders.begin(), end = builders.end();
	while(beg != end) {
		/* get the next operator and its subsuming mutants */
		const std::string & oprt = beg->first;
		MSGBuilder * builder = beg->second;
		SuMutantSet & smut = core->get_mutants(oprt);

		/* link the graph and delete builder */
		builder->link(); 
		builder->close(); 
		delete builder;

		/* update the subsuming mutants */
		smut.update_mutants();

		beg++;	/* to the next operator */
	}
	builders.clear();
}
void SuMutantExperimentDriver::derive_operator_II() {
	CError error(CErrorType::Runtime, 
		"SuMutantExperimentDriver::derive_operator_II", 
		"Invalid access: operators have not been designed");
	CErrorConsumer::consume(error);
	exit(CErrorType::Runtime);
}
void SuMutantExperimentDriver::derive_global_III() {
	CError error(CErrorType::Runtime,
		"SuMutantExperimentDriver::derive_global_III",
		"Invalid access: operators have not been designed");
	CErrorConsumer::consume(error);
	exit(CErrorType::Runtime);
}

void SuMutantExperimentOutput::write(const SuMutantExperimentCore & core) {
	/* write _summary_.txt */
	std::ofstream out1(dir->get_path() + "/_summary_.txt");
	this->write_summary(core, out1); out1.close();

	/* write for each operator.txt */
	const std::set <std::string> & keys = core.get_keys();
	auto beg = keys.begin(), end = keys.end();
	while (beg != end) {
		const std::string & key = *(beg++);
		const SuMutantSet & smut = core.get_subsuming_mutants(key);
		
		std::ofstream out2(dir->get_path() + "/" + key + ".txt");
		this->write_mutants(smut, out2); out2.close();
	}
}
void SuMutantExperimentOutput::write_summary(
	const SuMutantExperimentCore & core, std::ostream & out) {
	out << "operator\tequivalence\tsubsuming mutants\tsubsuming clusters\ttotal\n";

	const std::set <std::string> & keys = core.get_keys();
	auto beg = keys.begin(), end = keys.end();
	while (beg != end) {
		/* get next operator and its mutants */
		const std::string & oprt = *(beg++);
		const SuMutantSet & smut = core.get_subsuming_mutants(oprt);
		/* get the mutants by their categories */
		const MutantSet * eqs = smut.get_equivalent_mutants();
		const std::set<MuCluster *> & sclusters = smut.get_subsuming_clusters();
		const MutantSet & smutants = smut.get_subsuming_mutants();
		MSGraph & graph = smut.get_graph();

		/* output line */
		out << oprt << "\t";
		if (eqs == nullptr) out << 0;
		else out << eqs->number_of_mutants();
		out << "\t" << smutants.number_of_mutants() << "\t";
		out << sclusters.size() << "\t";
		out << graph.get_mutants().number_of_mutants() << "\n";
	}
}
void SuMutantExperimentOutput::write_mutants(
	const SuMutantSet & smut, std::ostream & out) {
	/* validate */
	if (dir == nullptr) {
		CError error(CErrorType::Runtime,
			"SuMutantExperimentOutput::write_mutants", "No file is opened");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}

	/* get the mutants in the graph */
	MSGraph & graph = smut.get_graph();
	MutantSpace & mspace = graph.get_space();
	const MutantSet & mutants = graph.get_mutants();
	const TextBuild & text = *(mspace.get_code_file().get_text());

	/* title */
	out << "mutant\tcluster\toperator\tline\torigin\treplace\tcategory\n";

	/* iterate each mutant in the space */
	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	while (mid < num) {
		/* get the next mutant in the graph */
		if (mutants.has_mutant(mid)) {
			/* get next mutant */
			Mutant & mutant = mspace.get_mutant(mid);
			MuCluster & cluster = graph.get_cluster_of(mid);
			const std::string & oprt = mutant.get_operator();
			const Mutation & mutation = 
				mutant.get_mutation(mutant.get_orders() - 1);
			const CodeLocation & loc = mutation.get_location();
			size_t line = text.lineOfIndex(loc.get_bias());
			std::string origin = loc.get_text_at();
			std::string replace = mutation.get_replacement();
			trim(origin); trim(replace);

			/* output line */ 
			out << mutant.get_id() << "\t" << cluster.get_id() << "\t";
			out << oprt << "\t" << line << "\t" << origin << "\t" << replace << "\t";

			/* category */
			switch (smut.get_category_of(mid)) {
			case SuMutantSet::equivalent:
				out << "equivalent";	break;
			case SuMutantSet::subsuming:
				out << "subsuming";		break;
			case SuMutantSet::subsumed:
				out << "subsumed";		break;
			default:
				CError error(CErrorType::Runtime, 
					"SuMutantExperimentOutput::write_mutants", "Unknown category.");
				CErrorConsumer::consume(error); exit(CErrorType::Runtime);
			}
			out << "\n";
		}

		/* to the next one */
		mid = mid + 1; 
	} /* end while */

	/* end */ out << std::endl;
}
void SuMutantExperimentOutput::trim(std::string & str) {
	std::string buffer;
	int k = 0, n = str.length();
	while (k < n) {
		char ch = str[k++];
		if (ch != '\n')
			buffer += ch;
	}
	str = buffer;
}


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

/* test main */
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "schedule";
	TestType ttype = TestType::schedule;

	// get root file and analysis dir 
	File & root = *(new File(prefix + prname));
	const std::vector<File *> & files = root.list_files();
	auto fbeg = files.begin(), fend = files.end(); 
	File * dir = nullptr;
	while (fbeg != fend) {
		File * file = *(fbeg++);
		if (file->get_local_name() == "analysis") {
			dir = file; break;
		}
	}
	if (dir == nullptr) {
		CError error(CErrorType::Runtime, "main()", "./analysis is undefined");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}

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

		// create experiment platform for this code file
		SuMutantExperimentCore core(mspace);
		SuMutantExperimentDriver driver;
		SuMutantExperimentOutput output;

		// get score vector producer | consumer
		ScoreSource & score_src = cscore.get_source(cfile);
		ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
		ScoreProducer producer(score_func); ScoreConsumer consumer(score_func);
		
		// execute the experiment 
		driver.start(core);
		driver.derive_operator_I(producer, consumer);
		driver.finish();

		// print the experiment result.
		output.open(*dir);
		output.write(core);
		output.close();

		std::cout << "Experiment finished: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &cscore;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	// exit
	std::cout << "Press any key to exit...\n"; getchar(); exit(0);
}
