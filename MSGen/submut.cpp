#include "submut.h"

void MutGroup::clear() {
	/* clear equivalents */
	equivalents = nullptr;
	/* clear subsuming mutants */
	subsumings.clear();
	/* re-open group and clear graph */
	builder.open(graph); 
}
void MutGroup::add_mutant(Mutant::ID mid, const BitSeq & score_vector) {
	builder.add(mid, score_vector);
}
void MutGroup::classify() {
	/* classify */ builder.link();

	/* declarations */
	const std::set<MuCluster *> & roots = graph.get_roots();
	auto beg = roots.begin(), end = roots.end(); MuCluster * cluster;
	
	/* get equivalents */
	equivalents = nullptr;
	if (roots.size() == 1) {
		cluster = *beg;
		if (cluster->get_score_degree() == 0)
			equivalents = cluster;
	}
	
	/* get subsuming clusters */
	subsumings.clear();
	if (equivalents == nullptr) {
		beg = roots.begin(), end = roots.end();
		while (beg != end) subsumings.insert(*(beg++));
	}
	else {
		const std::vector<MuSubsume> & edges 
			= equivalents->get_ou_port().get_edges();
		auto ebeg = edges.begin(), eend = edges.end();
		while (ebeg != eend) {
			const MuSubsume & edge = *(ebeg++);
			MuCluster & trg = edge.get_target();
			subsumings.insert(&trg);
		}
	}
	
	/* return */ return;
}
MutCategory MutGroup::category_of(Mutant::ID mid) const {
	if (graph.has_cluster_of(mid)) {
		MuCluster * cluster = 
			&(graph.get_cluster_of(mid));
		if (cluster == equivalents) 
			return Equvalent_Category;
		else if (subsumings.count(cluster) > 0) 
			return Subsuming_Category;
		else return Subsumed_Category;
	}
	else return NotBelong_Category;
}

const MutGroup & MutLevel::get_global_group() const {
	if (global_group == nullptr) {
		CError error(CErrorType::Runtime, 
			"MutLevel::get_global_group", 
			"Invalid access: group = nullptr");
		CErrorConsumer::consume(error); 
		exit(CErrorType::Runtime);
	}
	else return *global_group;
}
const MutGroup & MutLevel::get_operator_group(const std::string & op) const {
	if (operator_groups.count(op) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"MutLevel::get_operator_group", 
			"Invalid operator: " + op);
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = operator_groups.find(op);
		return *(iter->second);
	}
}
void MutLevel::add(Mutant::ID mid, const BitSeq & score_vector) {
	/* get operator name */
	Mutant & mutant = mspace.get_mutant(mid);
	const std::string oprt = mutant.get_operator();

	/* validate the operator name */
	if (valid_operator(oprt)) {
		/* get mutant operator group */
		MutGroup * group;
		if (operator_groups.count(oprt) == 0) {
			group = new MutGroup(mspace);
			operator_groups[oprt] = group;
			operators.insert(oprt);
		}
		else {
			auto iter = operator_groups.find(oprt);
			group = iter->second;
		}

		/* add mutant */ 
		group->add_mutant(mid, score_vector);
	}
}
void MutLevel::end() {
	/* extract subsuming mutants for each operator */
	auto beg = operator_groups.begin();
	auto end = operator_groups.end();
	while (beg != end) {
		MutGroup * group = (beg++)->second;
		group->classify();
	}

	/* extract subsuming mutants from all operators' subsuming mutants */
	global_group = new MutGroup(mspace);
	beg = operator_groups.begin();
	end = operator_groups.end();
	while (beg != end) {
		MutGroup * group = (beg++)->second;
		const std::set<MuCluster *> & 
			subsuming_clusters = group->get_subsumings();

		auto cbeg = subsuming_clusters.begin();
		auto cend = subsuming_clusters.end();
		while (cbeg != cend) {
			MuCluster * cluster = *(cbeg++);
			if (cluster->get_mutants().number_of_mutants() != 0) {
				Mutant::ID mid = get_mutant_of(*cluster);
				const BitSeq & svec = cluster->get_score_vector();
				global_group->add_mutant(mid, svec);
			}
		}
	}
	global_group->classify();

	/* return */ return;
}
void MutLevel::clear() {
	if(global_group != nullptr)
		delete global_group;
	global_group = nullptr;

	auto beg = operator_groups.begin();
	auto end = operator_groups.end();
	while (beg != end) 
		delete (beg++)->second;
	operator_groups.clear();
}
Mutant::ID MutLevel::get_mutant_of(const MuCluster & cluster) {
	const MutantSet & mutants = cluster.get_mutants();
	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	while (mid < num) {
		if (mutants.has_mutant(mid))
			return mid;
		else mid++;
	}

	CError error(CErrorType::Runtime, 
		"MutLevel::get_mutant_of", 
		"Empty cluster");
	CErrorConsumer::consume(error);
	exit(CErrorType::Runtime);
}
static std::set<std::string> debug_operators;
bool MutLevel::valid_operator(const std::string & op) {
	if (op == "I-CovAllNod" || op == "I-CovAllEdg" || op == "u-VDTR"
		|| op == "u-VTWD" || op == "I-DirVarIncDec" || op == "I-IndVarIncDec" || op == "u-Oido"
		|| op == "I-DirVarAriNeg" || op == "I-DirVarBitNeg" || op == "I-DirVarLogNeg"
		|| op == "I-IndVarAriNeg" || op == "I-IndVarBitNeg" || op == "I-IndVarLogNeg"
		|| op == "II-ArgAriNeg" || op == "II-ArgBitNeg" || op == "II-ArgLogNeg"
		|| op == "u-OCNG" || op == "u-OLNG"
		|| op == "u-SBRC" || op == "u-SCRB"
		|| op == "u-SRSR" || op == "u-SSDL"
		|| op == "u-SDWD" || op == "u-SWDD"
		/* new considered */
		|| op == "u-SMTC" || op == "u-SMTC" || op == "u-SMVB"
		|| op == "u-STRI" || op == "u-STRP" || op == "u-SSWM"
		/* end: new considered */
		|| op == "u-Cccr" || op == "u-Ccsr" || op == "u-CRCR"
		|| op == "u-VLSR" || op == "u-VGSR" || op == "u-VLPR" || op == "u-VGPR"
		|| op == "u-OAAN" || op == "u-OAAA" || op == "u-OABN" || op == "u-OABA" || op == "u-OAEA" || op == "u-OARN" || op == "u-OALN" || op == "u-OALA" || op == "u-OASN" || op == "u-OASA"
		|| op == "u-OBAN" || op == "u-OBAA" || op == "u-OBBN" || op == "u-OBBA" || op == "u-OBEA" || op == "u-OBRN" || op == "u-OBLN" || op == "u-OBLA" || op == "u-OBSN" || op == "u-OBSA"
		|| op == "u-OLAN" || op == "u-OLAA" || op == "u-OLBN" || op == "u-OLBA" || op == "u-OLEA" || op == "u-OLRN" || op == "u-OLLN" || op == "u-OLLA" || op == "u-OLSN" || op == "u-OLSA"
		|| op == "u-OSAN" || op == "u-OSAA" || op == "u-OSBN" || op == "u-OSBA" || op == "u-OSEA" || op == "u-OSRN" || op == "u-OSLN" || op == "u-OSLA" || op == "u-OSSN" || op == "u-OSSA"
		|| op == "u-ORAN" || op == "u-ORAA" || op == "u-ORBN" || op == "u-ORBA" || op == "u-OREA" || op == "u-ORRN" || op == "u-ORLN" || op == "u-ORLA" || op == "u-ORSN" || op == "u-ORSA"
		|| op == "u-OEAN" || op == "u-OEAA" || op == "u-OEBN" || op == "u-OEBA" || op == "u-OERN" || op == "u-OELN" || op == "u-OELA" || op == "u-OESN" || op == "u-OESA") return true;
	else {
		debug_operators.insert(op);
		return false;
	}
}

MutOutputer::MutOutputer(const File & root) {
	const std::vector<File *> & files = root.list_files();
	auto beg = files.begin(), end = files.end(); dir = nullptr;
	while(beg != end) {
		File & file = *(*(beg++));
		if (file.get_local_name() == "analysis") {
			dir = &file;
		}
	}

	if (dir == nullptr) {
		CError error(CErrorType::InvalidArguments, 
			"MutOutputer::MutOutputer", 
			"Undefined: ../analysis/");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
}
void MutOutputer::write(const MutLevel & data) {
	std::ofstream out1(dir->get_path() + "/summary.txt");
	write_summary(data, out1); out1.close();

	std::ofstream out2(dir->get_path() + "/mutants.txt");
	write_mutants(data, out2); out2.close();

	std::ofstream out3(dir->get_path() + "/distribution.txt");
	write_distribution(data, out3); out3.close();
}
void MutOutputer::write_summary(const MutLevel & data, std::ostream & out) {
	/* declaration */
	size_t M = 0, Eq = 0, MOP = 0, SOP = 0;
	size_t MOP_S = 0, SOP_S = 0, G_S = 0;

	MutantSpace & mspace = data.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	std::set<std::string> SOP_set; 
	const MutGroup & group = data.get_global_group();

	/* iterate the mutants by operators applied */
	const std::set<std::string> & OP_set = data.get_operators();
	auto beg = OP_set.begin(), end = OP_set.end();
	while (beg != end) {
		/* get next operator and its mutants group */
		const std::string & op = *(beg++);
		const MutGroup & op_group = data.get_operator_group(op);

		/* counts total summary */
		M += op_group.get_mutants().number_of_mutants();
		if (op_group.get_equivalents() != nullptr)
			Eq += ((op_group.get_equivalents())->size());

		/* count MOP and MOP_S */
		const std::set<MuCluster *> & clusters = op_group.get_subsumings();
		auto cbeg = clusters.begin(), cend = clusters.end();
		while (cbeg != cend) {
			MuCluster & cluster = *(*(cbeg++));
			MOP_S += cluster.size();
		}
		MOP = MOP + 1;
	}

	/* iterate global subsuming mutants for SOP_set */
	for (mid = 0; mid < num; mid++) {
		if (group.category_of(mid) == Subsuming_Category) {
			Mutant & mutant = mspace.get_mutant(mid);
			SOP_set.insert(mutant.get_operator()); G_S++;
		}
	}

	/* count SOP-subsuming mutants */
	beg = SOP_set.begin(), end = SOP_set.end();
	while (beg != end) {
		const std::string & op = *(beg++);
		const MutGroup & op_group = data.get_operator_group(op);

		/* count SOP-subsuming mutants */
		const std::set<MuCluster *> & clusters = op_group.get_subsumings();
		auto cbeg = clusters.begin(), cend = clusters.end();
		while (cbeg != cend) {
			MuCluster & cluster = *(*(cbeg++));
			SOP_S += cluster.size();
		}
		SOP = SOP + 1;
	}

	/* output */
	out << "Total mutants: \t" << M << "\n";
	out << "Equivalent set:\t" << Eq << "\n";
	out << "Total operator:\t" << MOP << "\n";
	out << "Subsuming ops: \t" << SOP << "\n";
	out << "MOP-Subsumings:\t" << MOP_S << "\n";
	out << "SOP-Subsumings:\t" << SOP_S << "\n";
	out << "Glb-Subsumings:\t" << G_S << "\n";

	out << "\nSubsuming Operators List:\n";
	beg = SOP_set.begin(), end = SOP_set.end();
	int count = 0;
	while (beg != end) {
		out << *(beg++) << "; ";
		if(count++ == 4) {
			count = 0; out << "\n";
		}
	}
	out << "\n";

	/* return */ out << std::endl; return;
}
void MutOutputer::write_mutants(const MutLevel & data, std::ostream & out) {
	/* declarations */
	MutantSpace & mspace = data.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	const MutGroup & global_group = data.get_global_group();
	const TextBuild & text = *(mspace.get_code_file().get_text());

	/* title */
	out << "mid\tline\torigin\treplace\toperator\tmode\top_cluster\top_category\tcluster\tcategory\n";

	/* output each mutant and their elements */
	for (mid = 0; mid < num; mid++) {
		/* get next mutant of valid operator */
		Mutant & mutant = mspace.get_mutant(mid);
		const std::string & oprt = mutant.get_operator();
		if (!data.has_operator(oprt)) continue;

		/* get the valid mutant in data space */
		const MutGroup & group = data.get_operator_group(oprt);
		if (!group.get_mutants().has_mutant(mid)) continue;
		
		/* basic information */
		const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & location = mutation.get_location();
		std::string origin = location.get_text_at();
		std::string replace = mutation.get_replacement();
		trim(origin); trim(replace);
		out << mutant.get_id() << "\t";
		out << text.lineOfIndex(location.get_bias()) << "\t";
		out << origin << "\t" << replace << "\t";

		/* operator mode */
		out << oprt << "\t" << "???" << "\t";

		/* operator category */
		MuCluster & op_cluster = group.get_cluster_of(mid);
		out << op_cluster.get_id() << "\t";
		switch (group.category_of(mid)) {
		case Equvalent_Category:
			out << "equivalent"; break;
		case Subsuming_Category:
			out << "subsuming"; break;
		case Subsumed_Category:
			out << "subsumed"; break;
		default:
			out << ""; break;
		}
		out << "\t";

		/* global test */
		if (global_group.get_mutants().has_mutant(mid)) {
			MuCluster & global_cluster 
				= global_group.get_cluster_of(mid);
			out << global_cluster.get_id() << "\t";

			switch (global_group.category_of(mid)) {
			case Equvalent_Category:
				out << "equivalent"; break;
			case Subsuming_Category:
				out << "subsuming"; break;
			case Subsumed_Category:
				out << "subsumed"; break;
			default:
				out << ""; break;
			}
		}

		out << "\n";	/* next line */
	} /* end for */

	/* end */
	out << std::endl; return;
}
void MutOutputer::write_distribution(const MutLevel & data, std::ostream & out) {
	/* count global-subsuming for each operator */
	const MutGroup & group = data.get_global_group();
	std::map<std::string, size_t> op_subsumings;
	MutantSpace & mspace = data.get_space(); size_t GS = 0;
	Mutant::ID mid, num = mspace.number_of_mutants();
	for (mid = 0; mid < num; mid++) {
		if (group.category_of(mid) == Subsuming_Category) {
			Mutant & mutant = mspace.get_mutant(mid);
			const std::string & op = mutant.get_operator();
			if (op_subsumings.count(op) == 0)
				op_subsumings[op] = 0;

			auto iter = op_subsumings.find(op);
			size_t count = iter->second; GS++;
			op_subsumings[op] = count + 1;
		}
	}

	/* print title */
	out << "Op\t#Mut\tEquiv\tOp-Su\tG-Su\tEfficiency\tContribution\n";
	const std::set<std::string> & operators = data.get_operators();
	
	/* print lines */
	auto obeg = operators.begin(), oend = operators.end();
	while (obeg != oend) {
		/* get next operator group */
		const std::string & op = *(obeg++);
		const MutGroup & op_group = data.get_operator_group(op);

		/* #Mut */
		size_t M = op_group.get_mutants().number_of_mutants();

		/* #Equiv */
		size_t E = 0;
		if (op_group.get_equivalents() != nullptr)
			E = op_group.get_equivalents()->size();

		/* #Op-Su */
		size_t OpS = op_group.get_subsumings().size();

		/* #G-Su */
		size_t GSu = 0;
		if (op_subsumings.count(op) > 0) {
			auto iter = op_subsumings.find(op);
			GSu = iter->second;
		}

		/* efficiency */
		double eff = ((double)GSu) / ((double)OpS);
		double ctr = ((double)GSu) / GS;

		/* output */
		out << op << "\t";
		out << M << "\t";
		out << E << "\t";
		out << OpS << "\t";
		out << GSu << "\t";
		out << eff << "\t";
		out << ctr << "\n";
	}
	op_subsumings.clear();

	/* return */ out << std::endl; return;
}
void MutOutputer::trim(std::string & str) {
	std::string bak; char ch;
	int i = 0, n = str.length();
	while (i < n) {
		ch = str[i++];

		if (ch == '\n') bak += " ";
		else if (ch == '\t') bak += " ";
		else bak += ch;
	}
	str = bak;
}

void TestMachine::start(const std::set<std::string> & ops) {
	close();
	auto beg = ops.begin(), end = ops.end();
	while (beg != end) {
		const std::string & op = *(beg++);
		if(context.has_operator(op))
			operators.insert(op);
	}
}
double TestMachine::evaluate(const TestSet & tests) {
	/* initialization */
	const MutGroup & group = context.get_global_group(); 
	const std::set<MuCluster *> & dom_mutants = group.get_subsumings();
	size_t K = 0, M = dom_mutants.size();

	/* validation */
	if (M == 0) {
		CError error(CErrorType::Runtime, "TestMachine::evaluate", 
			"Invalid context: empty global subsuming mutants");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}

	/* count the number of killed dominator mutants */
	auto beg = dom_mutants.begin(), end = dom_mutants.end();
	while (beg != end) {
		MuCluster & cluster = *(*(beg++));
		if (is_killed(tests, cluster)) {
			K++; 
		}
	}

	/* TODO evaluate information */
	std::cerr << "Dominator score: " << K << "/" << M << "\n";

	/* return */ return ((double) K) / ((double) M);
}
bool TestMachine::is_killed(const TestSet & tests, const MuCluster & cluster) {
	const BitSeq & tseq = tests.get_set_vector();
	const BitSeq & mseq = cluster.get_score_vector();
	BitSeq rseq(mseq); rseq.conjunct(tseq);
	//std::cerr << "\tResult: " << rseq.degree() << "\t" << rseq.all_zeros() << "\t" << !(rseq.all_zeros()) << "\n";
	return !(rseq.all_zeros());
}
void TestMachine::generate(TestSet & tests) {
	/* set requirements */
	std::set<MuCluster *> requirements;

	/* add all subsuming clusters into requirements */
	auto obeg = operators.begin();
	auto oend = operators.end();
	while (obeg != oend) {
		/* get the next operator-group */
		const std::string & oprt = *(obeg++);
		const MutGroup & op_group = 
			context.get_operator_group(oprt);
		const std::set<MuCluster *> &
			subsumings = op_group.get_subsumings();

		/* add the subsuming clusters */
		auto sbeg = subsumings.begin();
		auto send = subsumings.end();
		while (sbeg != send) 
			requirements.insert(*(sbeg++));
	}

	/* TODO output requirements Op-S: for debug */  
	std::cerr << "Select requirements: " << requirements.size() << "\n";

	/* generate template */
	std::vector<BitSeq *> templates; std::vector<unsigned> seeds;
	this->select_minimal_template(requirements, templates, tests.get_space());

	/* generate test seeds */
	for (int i = 0; i < templates.size(); i++) 
		seeds.push_back(1);

	/* generate test suite */
	this->generate_test_suite(templates, seeds, tests);

	/* clear test templates */
	for (int i = 0; i < templates.size(); i++)
		delete templates[i];
	seeds.clear(); templates.clear();

	/* return */ return;
}
void TestMachine::generate_test_suite(
	const std::vector<BitSeq *> & templates,
	const std::vector<unsigned> & seeds,
	TestSet & tests) {
	/* validation & initialization */
	if (seeds.size() != templates.size()) {
		CError error(CErrorType::InvalidArguments, 
			"TestMachine::generate_test_suite", 
			"Inconsistent inputs");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else tests.clear();

	/* insert test id into set */
	for (int i = 0; i < seeds.size(); i++) {
		unsigned int seed = seeds[i];
		const BitSeq & ri = *(templates[i]);
		TestCase::ID tid = this->find_test_at(ri, seed);
		tests.add_test(tid);
	}

	/* return */ return;
}
TestCase::ID TestMachine::find_test_at(const BitSeq & bits, TestCase::ID seed) {
	TestCase::ID k, n = bits.bit_number(); 
	unsigned origin_seed = seed;
	for (k = 0; k < n; k++) {
		if (bits.get_bit(k) == BIT_1) {
			if ((--seed) == 0)
				return k;
		}
	}

	/* not found */
	CError error(CErrorType::InvalidArguments, "TestMachine::find_test_at", 
		"Undefined index (" + std::to_string(origin_seed) + ")");
	CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
}
void TestMachine::select_minimal_template(
	const std::set<MuCluster *> & requirements,
	std::vector<BitSeq *> & templates,
	const TestSpace & test_space) {
	/* clusters to be eliminated */ 
	std::set<MuCluster *> ER, RS; templates.clear();
	BitSeq bits1(test_space.number_of_tests());
	BitSeq bits2(test_space.number_of_tests());

	/* initialization */
	auto beg = requirements.begin();
	auto end = requirements.end();
	while (beg != end) RS.insert(*(beg++));

	/* compute template for each requirement */
	while (!RS.empty()) {
		/* cluster the next groups */
		beg = RS.begin(), end = RS.end();
		
		/* initialize test suite */
		MuCluster * next = *(beg++);
		bits1.assign(next->get_score_vector());
		bits2.assign(next->get_score_vector());

		/* compute clusters with common tests (maximum) */
		ER.insert(next);
		while (beg != end) {
			next = *(beg++); 
			bits2.conjunct(next->get_score_vector());
			if (!bits2.all_zeros()) {
				ER.insert(next);
				bits1.assign(bits2);
			}
			else {
				bits2.assign(bits1);
			}
		}

		/* add new test set */
		BitSeq * seq = new BitSeq(bits1);
		templates.push_back(seq);

		/* eliminate ER */
		beg = ER.begin(), end = ER.end();
		while (beg != end) 
			RS.erase(*(beg++));
	}

	/* return */ return;
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
/* generate sufficient operators for traditional selective mutations */
static void selective_operators(std::set<std::string> & oprts) {
	oprts.clear();

	oprts.insert("u-SWDD");
	oprts.insert("u-SMTC");
	oprts.insert("u-SSDL");
	oprts.insert("u-OLBN");
	oprts.insert("u-OASN");
	oprts.insert("u-ORRN");
	oprts.insert("u-VTWD");
	oprts.insert("u-VDTR");
	oprts.insert("u-Cccr");
	oprts.insert("u-Ccsr");
}
/* generate operators based on subsuming mutations */
static void subsuming_operators(std::set<std::string> & oprts) {
	oprts.clear();

	oprts.insert("u-VDTR"); oprts.insert("u-VTWD");
	oprts.insert("u-VLSR"); oprts.insert("u-VGSR");
	oprts.insert("u-OAAN"); oprts.insert("u-OABN");
	oprts.insert("u-OARN"); oprts.insert("u-ORAN");
	oprts.insert("u-ORRN"); oprts.insert("u-ORSN");
	oprts.insert("I-DirVarIncDec");
	oprts.insert("I-IndVarIncDec");
}

/* test */

int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "schedule2"; 
	TestType ttype = TestType::schedule2;

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

		// get score vector producer | consumer
		ScoreSource & score_src = cscore.get_source(cfile);
		ScoreFunction & score_func = *(score_src.create_function(tests, mutants));
		ScoreProducer producer(score_func); ScoreConsumer consumer(score_func);

		// get experiment analyzer and data 
		MutAnalyzer analyzer(mspace);
		debug_operators.clear();
		analyzer.update(producer, consumer);
		MutLevel & data = analyzer.get_data();
		std::cout << "Compute subsuming mutants: finished\n";

		// end output
		MutOutputer outputer(root);
		outputer.write(data);
		std::cout << "Output subsuming mutants: finished\n";

		// evaluation 
		std::set<std::string> oprts;
		// selective_operators(oprts);
		subsuming_operators(oprts);

		MutantSet & smutants = *(mspace.create_set());
		TestSet & stests = *(ctest.malloc_test_set());
		TestMachine machine(data);
		machine.start(oprts);
		std::cout << "\nStart to generate test...\n";
		machine.generate(stests);
		std::cout << "Select " << stests.size() << " tests from set...\n";
		double score = machine.evaluate(stests);
		std::cout << "Score is " << score * 100 << "%...\n";

	}

	// delete memory
	delete &cscore;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	// debug operators
	auto cbeg = debug_operators.begin();
	auto cend = debug_operators.end();
	std::cout << "\nUnused operators:\n";
	while (cbeg != cend) {
		const std::string & op = *(cbeg++);
		std::cout << op << "\n";
	}
	std::cout << "\n";

	// exit
	std::cout << "Press any key to exit...\n"; getchar(); exit(0);
}
