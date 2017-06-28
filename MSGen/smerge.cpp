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

	/* compute the subsuming mutant (without duplicated) in parent graph */
	records.clear();
	while (mid < num) {
		if (mutants.has_mutant(mid)) {
			/* get the cluster of this mutant */
			MuCluster & cluster = graph.get_cluster_of(mid);

			/* to avoid duplicated subsuming mutant */
			if (records.count(&cluster) == 0) {
				/* add subsuming mutant (without duplication) to parent graph */
				const BitSeq & bits = cluster.get_score_vector();
				builder.add(mid, bits);

				/* record the cluster */ records.insert(&cluster);
			}
		}
		/* to the next mutant */ mid = mid + 1;
	}
	records.clear();
}
void SuMutantMerger::extract() {
	builder.link(); 
	builder.close();
	answer->update_mutants();
}

SuMutantSet & SuMutantExperimentCore::get_global_mutants() {
	if (global_mut == nullptr) {
		MSGraph * graph = new MSGraph(mspace);
		global_mut = new SuMutantSet(GLOBAL, *graph);
	}
	return *global_mut;
}
SuMutantSet & SuMutantExperimentCore::get_mutants_I(const std::string & oprt) {
	SuMutantSet * smut;
	if (mut_I.count(oprt) == 0) {
		MSGraph * graph = new MSGraph(mspace);
		smut = new SuMutantSet(oprt, *graph);
		mut_I[oprt] = smut;
	}
	else {
		auto iter = mut_I.find(oprt);
		smut = iter->second;
	}
	return *smut;
}
SuMutantSet & SuMutantExperimentCore::get_mutants_II(const std::string & oprt) {
	SuMutantSet * smut;
	if (mut_II.count(oprt) == 0) {
		MSGraph * graph = new MSGraph(mspace);
		smut = new SuMutantSet(oprt, *graph);
		mut_II[oprt] = smut;
	}
	else {
		auto iter = mut_II.find(oprt);
		smut = iter->second;
	}
	return *smut;
}
SuMutantSet & SuMutantExperimentCore::get_mutants_III(const std::string & oprt) {
	SuMutantSet * smut;
	if (mut_III.count(oprt) == 0) {
		MSGraph * graph = new MSGraph(mspace);
		smut = new SuMutantSet(oprt, *graph);
		mut_III[oprt] = smut;
	}
	else {
		auto iter = mut_III.find(oprt);
		smut = iter->second;
	}
	return *smut;
}
void SuMutantExperimentCore::clear_mutants() {
	/* clear the subsuming mutants in level-I */
	auto beg = mut_I.begin(), end = mut_I.end();
	while (beg != end) {
		SuMutantSet * smut = (beg++)->second;
		MSGraph & graph = smut->get_graph();
		delete smut;
		delete &graph;
	}

	/* clear the subsuming mutants in level-II */
	beg = mut_II.begin(), end = mut_II.end();
	while (beg != end) {
		SuMutantSet * smut = (beg++)->second;
		MSGraph & graph = smut->get_graph();
		delete smut;
		delete &graph;
	}

	/* clear the subsuming mutants in level-III */
	beg = mut_III.begin(), end = mut_III.end();
	while (beg != end) {
		SuMutantSet * smut = (beg++)->second;
		MSGraph & graph = smut->get_graph();
		delete smut;
		delete &graph;
	}

	/* clear global subsuming mutant set */
	if (global_mut != nullptr) {
		MSGraph & graph = global_mut->get_graph();
		delete global_mut;
		delete &graph;
		global_mut = nullptr;
	}
	
	/* return */ return;
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
		SuMutantSet & smut = core->get_mutants_I(oprt);
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
		SuMutantSet & smut = core->get_mutants_I(oprt);

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
	
	const std::map<std::string, SuMutantSet *> 
		& mutants_I = core->get_mutants_I();
	auto beg1 = mutants_I.begin(), end1 = mutants_I.end();
	while (beg1 != end1) {
		/* get the next operator and its subsuming set */
		const std::string & oprt = beg1->first;
		SuMutantSet * smut = beg1->second; beg1++;

		/* get the operator name for level-II */
		type_II(oprt, type2);
		if (type2.empty()) continue;
		SuMutantSet & smut2 = core->get_mutants_II(type2);

		/* get the the append list and merger */
		SuMutantMerger * merger;
		if (mergers.count(type2) == 0) {
			merger = new SuMutantMerger();
			mergers[type2] = merger;
			merger->open(smut2);
		}
		else {
			auto iter = mergers.find(type2);
			merger = iter->second;
		}

		/* add the subsuming set of operator to this */
		merger->append(*smut);
	}

	/* initialize type-II-cache for constructing */
	auto mbeg = mergers.begin();
	auto mend = mergers.end();
	while (mbeg != mend) {
		/* get the next merger and its append-list */
		const std::string & type2 = mbeg->first;
		SuMutantMerger * merger = mbeg->second;

		/* create subsuming set for the type-2 */
		merger->extract(); merger->close();

		/* delete merger resource */
		delete merger; mbeg++;	
	}

	/* return */ return;
}
void SuMutantExperimentDriver::derive_operator_III() {
	/* declarations */
	std::string type2, type3; SuMutantSet * smut;
	std::map<std::string, SuMutantMerger *> mergers;
	SuMutantMerger * merger;

	/* iterate subsuming mutants at level-II */
	const std::map <std::string, SuMutantSet *> 
		& mut_II = core->get_mutants_II();
	auto beg = mut_II.begin(), end = mut_II.end();
	while (beg != end) {
		/* get next subsuming mutants for operators at level-II */
		type2 = beg->first; smut = beg->second; beg++;

		/* get the type name at level-III */
		type_III(type2, type3);
		if (type3.empty()) continue;

		/* get the type-III mutants */
		SuMutantSet & mut3 = core->get_mutants_III(type3);

		/* get the merger for the subsuming mutants */
		merger = nullptr;
		if (mergers.count(type3) == 0) {
			merger = new SuMutantMerger();
			merger->open(mut3);
			mergers[type3] = merger;
		}
		else {
			auto iter = mergers.find(type3);
			merger = iter->second;
		}

		/* append the mutants at level-II to merger */
		merger->append(*smut);
	}

	/* merge all the mutants at level-III */
	auto mbeg = mergers.begin();
	auto mend = mergers.end();
	while (mbeg != mend) {
		type3 = mbeg->first;
		merger = mbeg->second;
		mbeg++;

		merger->extract(); 
		merger->close(); 
		delete merger;
	}
	mergers.clear();
}
void SuMutantExperimentDriver::derive_operator_IV() {
	/* declarations */
	const std::map<std::string, SuMutantSet *> mut_III = core->get_mutants_III();
	auto beg = mut_III.begin(), end = mut_III.end(); SuMutantMerger merger;

	/* initial global graph */
	merger.open(core->get_global_mutants());
	/* append the mutants in child sets */
	while (beg != end) {
		SuMutantSet * mut3 = (beg++)->second;
		merger.append(*mut3);
	}
	/* extract the subsuming mutants from children */
	merger.extract(); merger.close();

	/* return */ return;
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
	else if (oprt == "u-SRSR" || oprt == "u-SWDD" || oprt == "u-SCRB" || oprt == "u-SBRC" || oprt == "I-RetStaRep")
		type2 = REP_STATEMENT;
	/* operand-error */
	else if (oprt == "u-Ccsr" || oprt == "u-CRCR" || oprt == "u-VSCR")
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
		|| oprt == "u-OLEA" || oprt == "u-OLLA" || oprt == "u-OLLN" || oprt == "u-OLRN"
		|| oprt == "u-OLSA" || oprt == "u-OLSN")
		type2 = OLXX;
	else if (oprt == "u-OSAN" || oprt == "u-OSAA" || oprt == "u-OSBN" || oprt == "u-OSBA"
		|| oprt == "u-OSEA" || oprt == "u-OSLA" || oprt == "u-OSLN" || oprt == "u-OSRN"
		|| oprt == "u-OSSA" || oprt == "u-OSSN")
		type2 = OSXX;
	else if (oprt == "u-ORAN" || oprt == "u-ORBN" || oprt == "u-ORRN" || oprt == "u-ORSN" || oprt == "u-ORLN")
		type2 = ORXX;
	/* ignored operator */
	else if (oprt == "I-DirVarRepCon" || oprt == "I-DirVarRepReq" || oprt == "II-ArgRepReq"
		|| oprt == "II-FunCalDel" || oprt == "I-IndVarRepCon" || oprt == "I-IndVarRepReq"
		|| oprt == "I-DirVarRepExt" || oprt == "I-DirVarRepGlo" || oprt == "I-DirVarRepLoc"
		|| oprt == "I-DirVarRepPar" || oprt == "I-IndVarRepExt" || oprt == "I-IndVarRepGlo"
		|| oprt == "I-IndVarRepLoc" || oprt == "I-IndVarRepPar" || oprt == "u-SMTC"
		|| oprt == "u-SMTT" || oprt == "u-SMVB" || oprt == "u-SSWM" || oprt == "u-OCOR"
		|| oprt == "II-ArgDel" || oprt == "u-SBRn" || oprt == "u-SGLR")
		type2 = "";
	/* unknown list */
	else {
		CError error(CErrorType::InvalidArguments, 
			"SuMutantExperimentDriver::type_II", 
			"Unknown operator: \"" + oprt + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}

	/* return */ return;
}
void SuMutantExperimentDriver::type_III(const std::string & type2, std::string & type3) {
	if (type2 == TRAP_STATEMENT || type2 == TRAP_CONDITION || type2 == TRAP_VALUE)
		type3 = TRAP;
	else if (type2 == INCDEC_REFER || type2 == INCDEC_VALUE)
		type3 = INCDEC;
	else if (type2 == NEG_BOOLEAN || type2 == NEG_BINARYS || type2 == NEG_NUMBERS)
		type3 = NEG;
	else if (type2 == OAXX || type2 == OBXX || type2 == OEXX || type2 == OLXX || type2 == ORXX || type2 == OSXX)
		type3 = OPERATOR;
	else if (type2 == VAR_TO_CONST || type2 == CONST_TO_CONST || type2 == VAR_TO_VAR)
		type3 = OPERAND;
	else if (type2 == DEL_STATEMENT || type2 == REP_STATEMENT)
		type3 = STATEMENT;
	else {
		CError error(CErrorType::InvalidArguments,
			"SuMutantExperimentDriver::type_III",
			"Unknown operator: \"" + type2 + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
}

void SuMutantExperimentOutput::write(const SuMutantExperimentCore & core) {
	/* write prevalence.txt */
	std::ofstream out1(dir->get_path() + "/prevalence.txt");
	this->write_prevalence(core, out1); out1.close();

	/* write mutations.txt */
	std::ofstream out2(dir->get_path() + "/mutations.txt");
	this->write_mutations(core, out2); out2.close();
}
void SuMutantExperimentOutput::write_mutations(const SuMutantExperimentCore & core, std::ostream & out) {
	/* get mutant sets in core data */
	const std::map<std::string, SuMutantSet *> & mut_I = core.get_mutants_I();
	const std::map<std::string, SuMutantSet *> & mut_II = core.get_mutants_II();
	const std::map<std::string, SuMutantSet *> & mut_III = core.get_mutants_III();
	const SuMutantSet & global_mut = core.get_mutants_IV();

	/* for iterating mutants in the space */
	MutantSpace & mspace = core.get_space();
	Mutant::ID k = 0, mid, num = mspace.number_of_mutants();
	const TextBuild & text = *(mspace.get_code_file().get_text());
	std::string origin, replace, type1, type2, type3;

	/* title for table */
	out << "mid\tline\torigin\treplace\ttyI\tcid\tcat\ttyII\tcid\tcat\ttyIII\tcid\tcat\tcid\tcat\n";

	/* iterate each mutant in the space */
	while (k < num) {
		/* get the next mutant in the space */
		Mutant & mutant = mspace.get_mutant(k++);
		mid = mutant.get_id();

		/* it's one of the mutant in some type-I cluster */
		type1 = mutant.get_operator();
		if (mut_I.count(type1) > 0) {
			/* get the subsuming set where the mutant is defined */
			auto iter1 = mut_I.find(type1);
			SuMutantSet & smut1 = *(iter1->second);

			if (smut1.has_mutant(mid)) {
				/* get its content */
				const Mutation & mutation =
					mutant.get_mutation(mutant.get_orders() - 1);
				const CodeLocation & loc = mutation.get_location();
				origin = loc.get_text_at();
				replace = mutation.get_replacement();
				trim(origin); trim(replace);

				/* base head for mutant */
				out << mutant.get_id() << "\t";
				out << text.lineOfIndex(loc.get_bias()) << "\t";
				out << origin << "\t" << replace << "\t";

				/* print information of category in level-I */
				out << type1 << "\t";
				out << smut1.get_cluster_of(mid).get_id() << "\t";
				out << get_category(smut1, mid) << "\t";

				/* to parent level */
				SuMutantExperimentDriver::type_II(type1, type2);
				if (!type2.empty() && mut_II.count(type2)) {
					/* get the subsuming set of level-II where the mutant is defined */
					auto iter2 = mut_II.find(type2);
					SuMutantSet & smut2 = *(iter2->second);

					if (smut2.has_mutant(mid)) {
						/* print information of category in level-II */
						out << type2 << "\t";
						out << smut2.get_cluster_of(mid).get_id() << "\t";
						out << get_category(smut2, mid) << "\t";

						/* to parent level */
						SuMutantExperimentDriver::type_III(type2, type3);
						if (!type3.empty() && mut_III.count(type3)) {
							/* get the subsuming set of level-III */
							auto iter3 = mut_III.find(type3);
							SuMutantSet & smut3 = *(iter3->second);

							if (smut3.has_mutant(mid)) {
								/* print information of category in level-III */
								out << type3 << "\t";
								out << smut3.get_cluster_of(mid).get_id() << "\t";
								out << get_category(smut3, mid) << "\t";

								/* to parent set */
								if (global_mut.has_mutant(mid)) {
									out << global_mut.get_cluster_of(mid).get_id() << "\t";
									out << get_category(global_mut, mid);
								}	/* end if: global */
							}
						}
					}
				}

				/* to the next line */ out << "\n";
			}
		}
	} /* end while */

	/* #EOF */ out << std::endl; 
}
void SuMutantExperimentOutput::write_prevalence(const SuMutantExperimentCore & core, std::ostream & out) {
	// do nothing...
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
std::string SuMutantExperimentOutput::get_category(const SuMutantSet & smut, Mutant::ID mid) {
	switch (smut.get_category_of(mid)) {
	case SuMutantSet::equivalent:	return "equivalent";
	case SuMutantSet::subsuming:	return "subsuming";
	case SuMutantSet::subsumed:		return "subsumed";
	case SuMutantSet::not_belong:	return "";
	default:
		CError error(CErrorType::Runtime, 
			"SuMutantExperimentOutput::get_category", 
			"Unknown category");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}
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
/* global subsuming mutants */
static void compute_subsumings(const CodeFile & cfile, MutantSet & mutants, TestSet & tests, CScore & cscore) {
	// get score vector producer | consumer
	ScoreSource & score_src = cscore.get_source(cfile);
	ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
	ScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

	MSGraph graph(mutants.get_space());
	MSGBuilder builder; ScoreVector * score_vector;

	builder.open(graph);
	while ((score_vector = producer.produce()) != nullptr) {
		Mutant::ID mid = score_vector->get_mutant();
		const BitSeq & bits = score_vector->get_vector();
		builder.add(mid, bits);
		consumer.consume(score_vector);
	}
	builder.link(); 
	builder.close();

	SuMutantSet smut("global", graph);
	std::cerr << "Subsuming Mutants: " << smut.get_subsuming_mutants().number_of_mutants()
		<< " {" << smut.get_subsuming_clusters().size() << "}\n";
}

/* test main */
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "bubble";
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
		driver.derive_operator_III();
		driver.derive_operator_IV();
		driver.finish();

		// print the experiment result.
		output.open(*dir);
		output.write(core);
		output.close();

		// test 
		//compute_subsumings(cfile, mutants, tests, cscore);
		std::cerr << "Subsuming Mutants: " << core.get_mutants_IV().get_subsuming_mutants().number_of_mutants()
			<< "{" << core.get_mutants_IV().get_subsuming_clusters().size() << "}\n";
		std::cout << "Experiment finished: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &cscore;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	// exit
	std::cout << "Press any key to exit...\n"; getchar(); exit(0);
}
