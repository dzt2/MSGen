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

MuCluster * SuOperatorSearch::belong_to(
	MuCluster * src, const std::set<MuCluster *> & SC) {
	const BitSeq & x = src->get_score_vector();
	auto beg = SC.begin(), end = SC.end();
	while (beg != end) {
		MuCluster & target = *(*(beg++));
		const BitSeq & y = target.get_score_vector();
		if (x.equals(y)) return &target;
	}
	return nullptr;
}

unsigned SuOperatorSearch::eval(const std::string & op) {
	return eval_mutants(op);
//	return eval_operats(op);
//	return eval_smutant(op);
}
unsigned SuOperatorSearch::eval_mutants(const std::string & op) {
	const MutGroup & group = context.get_operator_group(op);
	return group.get_mutants().number_of_mutants();
}
unsigned SuOperatorSearch::eval_operats(const std::string & op) {
	return 1;
}
unsigned SuOperatorSearch::eval_smutant(const std::string & op) {
	const MutGroup & group = context.get_operator_group(op);
	return group.get_subsumings().size();
}
void SuOperatorSearch::start() {
	finish();	/* clear state space */

	const MutGroup & group = context.get_global_group();
	const std::set<MuCluster *> & scs = group.get_subsumings();

	/* get two maps */
	const std::set<std::string> & ops = context.get_operators();
	auto obeg = ops.begin(), oend = ops.end();
	while (obeg != oend) {
		/* get the op-subsuming clusters */
		const std::string & op = *(obeg++);
		const MutGroup & op_group = context.get_operator_group(op);
		const std::set<MuCluster *> & op_sm = op_group.get_subsumings();

		/* build the items in two maps for this operator */
		auto sbeg = op_sm.begin(), send = op_sm.end();
		while (sbeg != send) {
			/* get the subsuming cluster for the op-cluster */
			MuCluster * src = *(sbeg++);
			MuCluster * trg = belong_to(src, scs);

			/* build maps */
			if (trg != nullptr) {
				/* build item space */
				if (op_scs.count(op) == 0)
					op_scs[op] = new std::set<MuCluster *>();
				if (sc_ops.count(trg) == 0)
					sc_ops[trg] = new std::set<std::string>();

				/* update items */
				auto iter1 = op_scs.find(op);
				(iter1->second)->insert(trg);
				auto iter2 = sc_ops.find(trg);
				(iter2->second)->insert(op);
			}
		}	/* end two maps for operator */

	}	/* end two maps for all operators */

	/* validation */
	auto sbeg = scs.begin(), send = scs.end();
	while (sbeg != send) {
		MuCluster * trg = *(sbeg++);
		if (sc_ops.count(trg) == 0) {
			CError error(CErrorType::Runtime, 
				"SuOperatorSearch::start", 
				"Incorrect context");
			CErrorConsumer::consume(error);
			exit(CErrorType::Runtime);
		}
	}	/* end validation */
}
void SuOperatorSearch::finish() {
	while (!_stack.empty()) 
		pop_stat();

	auto obeg = op_scs.begin();
	auto oend = op_scs.end();
	while (obeg != oend) 
		delete (obeg++)->second;
	op_scs.clear();

	auto sbeg = sc_ops.begin();
	auto send = sc_ops.end();
	while (sbeg != send)
		delete (sbeg++)->second;
	sc_ops.clear();
}
void SuOperatorSearch::solve(std::set<std::string> & SOP) {
	/* initialize */ 
	int_stat(); SOP.clear();

	unsigned minvalue = context.get_space().number_of_mutants();
	while (!_stack.empty()) {
		_Item & top = *(_stack.top());	/* get the top item state */

		/* valiate constraint of min-value */
		if (top.value >= minvalue) {
			pop_stat();
		}
		/* solution found and update SOP */
		else if (top.sc_state.empty()) {
			gen_solution(SOP);
			minvalue = top.value;

			// TODO debug
			std::cerr << "Found(" << minvalue << "): \t{";
			std::cerr << "SOP = " << SOP.size() << "; \t";
			std::cerr << "STC = " << _stack.size() << "; ";
			std::cerr << "}\n";

			pop_stat();
		}
		/* iteration finished at this selected item */
		else if (top.cursor >= top.op_state.size()) {	
			pop_stat();
		}
		/* continue to the next item */
		else {
			psh_stat(); top.cursor = top.cursor + 1;
		}
	} /* end while: stack solving */

	/* return */ return;
}
void SuOperatorSearch::int_stat() {
	/* clear and initialize state */
	while (!_stack.empty()) pop_stat();
	_Item * item = new _Item();
	

	/* for operator-candidates */
	auto obeg = op_scs.begin();
	auto oend = op_scs.end();
	while (obeg != oend) {
		const std::string & op = (obeg++)->first;
		item->op_state.push_back(op);
	}
	/* for subsuming mutants state */
	auto sbeg = sc_ops.begin();
	auto send = sc_ops.end();
	while (sbeg != send) {
		MuCluster * trg = (sbeg++)->first;
		item->sc_state.insert(trg);
	}
	/* initialize cursor to the next level */
	item->cursor = 0; item->value = 0;

	/* push & return */ _stack.push(item);
}
void SuOperatorSearch::pop_stat() {
	if (_stack.empty()) {
		CError error(CErrorType::Runtime, 
			"SuOperatorSearch::pop_stat", 
			"Incorrect access: empty state");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else {
		delete _stack.top();
		_stack.pop();
	}
}
void SuOperatorSearch::psh_stat() {
	if (_stack.empty()) {
		CError error(CErrorType::Runtime,
			"SuOperatorSearch::psh_stat",
			"Incorrect access: empty state");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else {
		/* get next operator to be selected */
		_Item & item = *(_stack.top());
		const std::string & op 
			= item.op_state[item.cursor];

		/* create new item */
		_Item & next = *(new _Item());

		/* compute subsuming clusters state */
		const std::set<MuCluster *> & origin = item.sc_state;
		auto op_iter = op_scs.find(op);
		std::set<MuCluster *> & reduce = *(op_iter->second);
		auto org_beg = origin.begin(), org_end = origin.end();
		while (org_beg != org_end) {
			MuCluster * trg = *(org_beg++);
			if (reduce.count(trg) == 0)
				next.sc_state.insert(trg);
		}

		/* compute the sufficient operators */
		std::set<std::string> sufficientOP;
		auto sub_beg = next.sc_state.begin();
		auto sub_end = next.sc_state.end();
		while (sub_beg != sub_end) {
			auto iter = sc_ops.find(*(sub_beg++));
			const std::set<std::string> & ops = *(iter->second);
			auto beg = ops.begin(), end = ops.end();
			while (beg != end) 
				sufficientOP.insert(*(beg++));
		}

		/* compute the operator candidates */
		for (int i = item.cursor + 1; 
			i < item.op_state.size(); i++) {
			const std::string & cop = item.op_state[i];
			if (sufficientOP.count(cop) > 0)
				next.op_state.push_back(cop);
		}

		/* compute value */
		next.value = item.value + eval(op);

		/* initialize cursor */ next.cursor = 0;

		/* push at last */ _stack.push(&next);
	}
}
void SuOperatorSearch::gen_solution(std::set<std::string> & SOP) {
	/* declarations */
	SOP.clear(); std::stack<_Item *> items;

	/* compute the operators */
	while (!_stack.empty()) {
		_Item & item = *(_stack.top());
		int key = item.cursor - 1;
		if (key >= 0) 
			SOP.insert(item.op_state[key]);
		_stack.pop(); items.push(&item);
	}

	/* recover the stack */
	while (!items.empty()) {
		_Item * item = items.top();
		items.pop(); _stack.push(item);
	}
}

void OperatorWriter::open(const File & root) {
	close();	/* close the original stream */

	/* get the ../analysis directory */
	this->dir = nullptr;
	const std::vector<File *> & files = root.list_files();
	for (int i = 0; i < files.size(); i++) {
		File & file = *(files[i]);
		if (file.get_local_name() == "analysis") {
			dir = &file; break;
		}
	}

	/* analysis dir is not found */
	if (dir == nullptr) {
		CError error(CErrorType::InvalidArguments, 
			"OperatorWriter::open", 
			"path ../analysis is not found");
		CErrorConsumer::consume(error); 
		exit(CErrorType::InvalidArguments);
	}
}
void OperatorWriter::close() {
	dir = nullptr; 
}
void OperatorWriter::write(MutLevel & data, const std::set<std::string> & SOP) {
	if (dir == nullptr) {
		CError error(CErrorType::Runtime,
			"OperatorWriter::write", 
			"writer is not opened");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	
	std::ofstream out1(dir->get_path() + "/summary.txt");
	this->write_summary(data, SOP, out1); out1.close();
	std::ofstream out2(dir->get_path() + "/contribution.txt");
	this->write_contribution(data, out2); out2.close();

	/* return */ return;
}
void OperatorWriter::write_summary(
	MutLevel & data, const std::set<std::string> & SOP, std::ostream & out) {
	/* declarations */
	size_t M = 0, E = 0, S = 0, O = 0, SO = 0, SM = 0, SE = 0, SSM = 0;

	/* operators applied */
	const std::set<std::string> & oprts = data.get_operators();
	auto beg = oprts.begin(), end = oprts.end();
	while (beg != end) {
		const std::string & oprt = *(beg++);
		const MutGroup & op_group = data.get_operator_group(oprt);
		M += op_group.get_mutants().number_of_mutants();
		if (op_group.get_equivalents() != nullptr)
			E += op_group.get_equivalents()->size();
	}

	/* subsuming set */
	const MutGroup & group = data.get_global_group();
	S = group.get_subsumings().size();

	/* operator-subsuming operator */
	O = oprts.size(); SO = SOP.size();

	/* mutants in subsuming operators */
	auto sbeg = SOP.begin(), send = SOP.end();
	while (sbeg != send) {
		const MutGroup & op_group = data.get_operator_group(*(sbeg++));
		SM += op_group.get_mutants().number_of_mutants();
		if (op_group.get_equivalents() != nullptr)
			SE += op_group.get_equivalents()->size();
		SSM += op_group.get_subsumings().size();
	}

	/* numbers output */
	out << "Mutants: \t" << M << "\n";
	out << "Equival: \t" << E << "\n";
	out << "SuMutat: \t" << S << "\n";
	out << "Operatr: \t" << O << "\n";
	out << "SubOprt: \t" << SO << "\n";
	out << "SOpMutt: \t" << SM << "\n";
	out << "SOpEquv: \t" << SE << "\n";
	out << "SOpSubM: \t" << SSM << "\n";

	/* subsuming operator print */
	int count = 0; out << "\n";
	sbeg = SOP.begin(), send = SOP.end();
	while (sbeg != send) {
		const std::string & op = *(sbeg++);
		out << op << "; ";
		if (++count == 5) {
			out << "\n"; count = 0;
		}
	}
	out << "\n";

	/* return */ out << std::endl; return;
}
void OperatorWriter::write_contribution(
	MutLevel & data, std::ostream & out) {
	/* declarations */
	const MutGroup & group = data.get_global_group();
	const std::set<MuCluster *> & Su = group.get_subsumings();
	const std::set<std::string> & ops = data.get_operators();
	size_t SuNum = Su.size();

	/* title */
	out << "Operator\tMutants\tEquiv\tNumber\tCop\tEop\n";

	auto beg = ops.begin(), end = ops.end();
	while (beg != end) {
		/* get next operator */
		const std::string & op = *(beg++);
		const MutGroup & op_group = data.get_operator_group(op);
		const std::set<MuCluster *> & opSu = op_group.get_subsumings();

		/* Mut, Eq, Ct, Ef, Num */
		size_t Mut = op_group.get_mutants().number_of_mutants();
		size_t Eqv = 0;
		if (op_group.get_equivalents() != nullptr)
			Eqv = op_group.get_equivalents()->size();

		/* count the contribution */
		size_t Num = 0;
		auto obeg = opSu.begin(), oend = opSu.end();
		while (obeg != oend) {
			const MuCluster & cluster = *(*(obeg++));
			if (belong_to(cluster, Su) != nullptr)
				Num++;
		}

		/* Ct, Ef */
		double Ct, Ef = 0.0;
		Ct = ((double)Num) / ((double)SuNum);
		if(Mut - Eqv != 0)
			Ef = ((double)Num) / ((double)(Mut - Eqv));

		/* print line */
		out << op << "\t";
		out << Mut << "\t";
		out << Eqv << "\t";
		out << Num << "\t";
		out << Ct << "\t";
		out << Ef << "\n";
	}

	out << std::endl; return;
}
MuCluster * OperatorWriter::belong_to(const MuCluster & cluster, const std::set<MuCluster *> & ops) {
	const BitSeq & x = cluster.get_score_vector();
	auto beg = ops.begin(), end = ops.end();
	while (beg != end) {
		MuCluster & target = *(*(beg++));
		const BitSeq & y = target.get_score_vector();
		if (x.equals(y)) return &target;
	}
	return nullptr;
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

		// solve subsuming operators
		std::set<std::string> SOP;
		SuOperatorSearch search(data);
		search.start(); search.solve(SOP); search.finish();
		std::cout << "Solving subsuming operators: " << SOP.size() << "\n";

		// end output
		OperatorWriter writer;
		writer.open(root);
		writer.write(data, SOP);
		writer.close();
		std::cout << "Output subsuming mutants: finished\n";

		// evaluation 
//		std::set<std::string> oprts;
//		selective_operators(oprts);
//		subsuming_operators(oprts);
//
//		MutantSet & smutants = *(mspace.create_set());
//		TestSet & stests = *(ctest.malloc_test_set());
//		TestMachine machine(data);
//		machine.start(oprts);
//		std::cout << "\nStart to generate test...\n";
//		machine.generate(stests);
//		std::cout << "Select " << stests.size() << " tests from set...\n";
//		double score = machine.evaluate(stests);
//		std::cout << "Score is " << score * 100 << "%...\n";

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
