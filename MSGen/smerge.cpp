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
	MutantSpace & mspace = mutants.get_space();
	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	
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
	/* declarations */
	std::string type2;
	std::map<std::string, SuMutantMerger *> mergers;
	std::map<std::string, std::list<SuMutantSet *> *> append_lists;
	
	const std::set<std::string> & keys = core->get_keys();
	auto kbeg = keys.begin(), kend = keys.end();
	while (kbeg != kend) {
		/* get the next operator and its subsuming set */
		const std::string & oprt = *(kbeg++);
		this->type_II(oprt, type2);
		if (type2.empty()) continue;

		/* get the the append list and merger */
		SuMutantMerger * merger;
		std::list<SuMutantSet *> * alist;
		if (mergers.count(type2) == 0) {
			merger = new SuMutantMerger();
			mergers[type2] = merger;
			alist = new std::list<SuMutantSet *>();
			append_lists[type2] = alist;
		}
		else {
			auto iter = mergers.find(type2);
			merger = iter->second;
			auto iter2 = append_lists.find(type2);
			alist = iter2->second;
		}

		/* add the subsuming set of operator to this */
		SuMutantSet & smut = core->get_mutants(oprt);
		alist->push_back(&smut);
	}

	/* initialize type-II-cache for constructing */
	type_II_cache.clear();
	auto mbeg = mergers.begin();
	auto mend = mergers.end();

	while (mbeg != mend) {
		/* get the next merger and its append-list */
		const std::string & type2 = mbeg->first;
		SuMutantMerger * merger = mbeg->second;
		auto iter = append_lists.find(type2);
		std::list<SuMutantSet *> * alist = iter->second;

		/* create subsuming set for the type-2 */
		SuMutantSet & smut = core->get_mutants(type2);
		
		/* add child-subsuming-mutants to this one */
		auto beg = alist->begin();
		auto end = alist->end();
		merger->open(smut);
		while (beg != end) {
			SuMutantSet * smut2 = *(beg++);
			merger->append(*smut2);
		}
		merger->extract(); merger->close();

		/* insert the type-2 to cache */
		type_II_cache.insert(type2);

		/* delete merger resource */
		delete merger; alist->clear(); delete alist;

		mbeg++;	/* to the next group */
	}

	/* return */ return;
}
void SuMutantExperimentDriver::derive_global_III() {
	CError error(CErrorType::Runtime,
		"SuMutantExperimentDriver::derive_global_III",
		"Invalid access: operators have not been designed");
	CErrorConsumer::consume(error);
	exit(CErrorType::Runtime);
}
void SuMutantExperimentDriver::type_II(const std::string & oprt, std::string & type2) {
	/* initial */ type2 = "";

	/* trap-error */
	if (oprt == "u-STRP" || oprt == "I-CovAllNod")
		type2 = TRAP_STATEMENT;
	else if (oprt == "u-STRI" || oprt == "I-CovAllEdg")
		type2 = TRAP_CONDITION;
	else if (oprt == "u-VDTR")
		type2 = TRAP_VALUE;
	/* inc-dec-error */
	else if (oprt == "u-VTWD" || oprt == "II-ArgIncDec")
		type2 = INCDEC_VALUE;
	else if (oprt == "u-Oido" || oprt == "I-DirVarIncDec" || oprt == "I-IndVarIncDec")
		type2 = INCDEC_REFER;
	/* negate-error */
	else if (oprt == "I-DirVarLogNeg" || oprt == "II-ArgLogNeg" || oprt == "I-IndVarLogNeg" || oprt == "u-OCNG" || oprt == "u-OLNG")
		type2 = NEG_BOOLEAN;
	else if (oprt == "I-DirVarBitNeg" || oprt == "II-ArgBitNeg" || oprt == "I-IndVarBitNeg" || oprt == "u-OBNG")
		type2 = NEG_BINARYS;
	else if (oprt == "I-DirVarAriNeg" || oprt == "II-ArgAriNeg" || oprt == "I-IndVarAriNeg" || oprt == "u-OANG")
		type2 = NEG_NUMBERS;
	/* statement-error */
	else if (oprt == "u-SSDL" || oprt == "I-RetStaDel")
		type2 = DEL_STATEMENT;
	else if (oprt == "u-SRSR" || oprt == "u-SWDD" || oprt == "u-SCRB" || oprt == "I-RetStaRep")
		type2 = REP_STATEMENT;
	/* operand-error */
	else if (oprt == "u-Ccsr" || oprt == "u-CRCR")
		type2 = VAR_TO_CONST;
	else if (oprt == "u-Cccr")
		type2 = CONST_TO_CONST;
	else if (oprt == "u-VGSR" || oprt == "u-VLPR" || oprt == "u-VLSR" || oprt == "II-ArgStcAli")
		type2 = VAR_TO_VAR;
	/* operator-error */
	else if (oprt == "u-OAAN" || oprt == "u-OAAA" || oprt == "u-OABN" || oprt == "u-OABA"
		|| oprt == "u-OAEA" || oprt == "u-OALA" || oprt == "u-OALN" || oprt == "u-OASN"
		|| oprt == "u-OASA" || oprt == "u-OARN")
		type2 = OAXX;
	else if (oprt == "u-OBAN" || oprt == "u-OBAA" || oprt == "u-OBBN" || oprt == "u-OBBA"
		|| oprt == "u-OBEA" || oprt == "u-OBLA" || oprt == "u-OBLN" || oprt == "u-OBSN"
		|| oprt == "u-OBSA" || oprt == "u-OBRN")
		type2 = OBXX;
	else if (oprt == "u-OEAA" || oprt == "u-OEBA" || oprt == "u-OESA" || oprt == "u-OELA")
		type2 = OEXX;
	else if (oprt == "u-OLAN" || oprt == "u-OLAA" || oprt == "u-OLBA" || oprt == "u-OLBN"
		|| oprt == "u-OLLA" || oprt == "u-OLLN" || oprt == "u-OLRN" || oprt == "u-OLSA" || oprt == "u-OLSN")
		type2 = OLXX;
	else if (oprt == "u-ORAN" || oprt == "u-ORBN" || oprt == "u-ORRN" || oprt == "u-ORSN" || oprt == "u-ORLN")
		type2 = ORXX;
	/* ignored operator */
	else if (oprt == "I-DirVarRepCon" || oprt == "I-DirVarRepReq" || oprt == "II-ArgRepReq"
		|| oprt == "II-FunCalDel" || oprt == "I-IndVarRepCon" || oprt == "I-IndVarRepReq"
		|| oprt == "I-DirVarRepExt" || oprt == "I-DirVarRepClo" || oprt == "I-DirVarRepLoc"
		|| oprt == "I-DirVarRepPar" || oprt == "I-IndVarRepExt" || oprt == "I-IndVarRepGlo"
		|| oprt == "I-IndVarRepLoc" || oprt == "I-IndVarRepPar" || oprt == "u-SMTC"
		|| oprt == "u-SMTT" || oprt == "u-SMVB" || oprt == "u-SSWM" || oprt == "u-OCOR"
		|| oprt == "II-ArgDel")
		type2 = "";
	/* unknown list */
	else {
		CError error(CErrorType::InvalidArguments, 
			"SuMutantExperimentDriver::type_II", 
			"Unknown operator: \"" + oprt + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
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
	std::string prname = "mid";
	TestType ttype = TestType::general;

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
		driver.derive_operator_II();
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
