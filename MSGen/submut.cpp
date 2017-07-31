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
		//debug_operators.insert(op);
		//return false;
		return true;
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

	  /* validate solution */
	if (SOP.empty()) {
		CError error(CErrorType::Runtime,
			"SuOperatorSearch::solve",
			"Solution is not found");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else return;
}
void SuOperatorSearch::re_solve(std::set<std::string> & SOP) {
	/* initialize */
	int_stat(); SOP.clear(); 
	size_t constrain = sc_ops.size();	/* number of remaining clusters */
	unsigned maxvalue = 0;	/* how many mutants are reduced */

	while (!_stack.empty()) {
		_Item & top = *(_stack.top());	/* get the top item state */

		/* validate constrain */
		if (top.sc_state.size() < constrain) {
			pop_stat();
		}
		/* validate a solution */
		else if (top.cursor >= top.op_state.size()) {
			if (top.value >= maxvalue) {
				maxvalue = top.value;
				gen_re_solution(SOP);
			}
			pop_stat();
		}
		/* further search child space */
		else {
			rpsh_stat(); top.cursor = top.cursor + 1;
		}
	} /* end while stack */
	
	/* validate solution */
	if (SOP.empty()) {
		CError error(CErrorType::Runtime, 
			"SuOperatorSearch::re_solve", 
			"Solution is not found");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else return;

}
void SuOperatorSearch::gr_solve(std::set<std::string> & SOP) {
	/* initialization */
	int_stat(); SOP.clear();
	_Item & root = *(_stack.top());
	root.op_state.clear();

	/* update root item */
	MuCluster * first = *(root.sc_state.begin());
	auto first_iter = sc_ops.find(first);
	std::set<std::string> & first_ops = *(first_iter->second);
	auto first_beg = first_ops.begin();
	auto first_end = first_ops.end();
	while (first_beg != first_end) 
		root.op_state.push_back(*(first_beg++));

	/* compute the best solution */
	unsigned minvalue = context.get_space().number_of_mutants();
	while (!_stack.empty()) {
		_Item & top = *(_stack.top());	/* get the top item state */

		if (top.value >= minvalue) {	/* violate limit */
			this->pop_stat();
		}
		else if (top.sc_state.empty()) {	/* find feasible solution */
			minvalue = top.value;
			this->pop_stat();
			this->gen_solution(SOP);
		}
		else if (top.cursor >= top.op_state.size()) {	/* out of range */
			this->pop_stat();
		}
		else {	/* construction solution further */
			gpsh_stat(); top.cursor = top.cursor + 1;
		}

	} /* end while solution */

	/* validate solution */
	if (SOP.empty()) {
		CError error(CErrorType::Runtime,
			"SuOperatorSearch::gr_solve",
			"Solution is not found");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else return;
}
void SuOperatorSearch::suf_solve(const std::set<std::string> & init,
	std::set<std::string> & SufSOP, double alpha) {
	/* initialization */
	SufSOP.clear(); init_stat(init);
	size_t limit = sc_ops.size() * alpha;
	unsigned maxvalue = 0;

	while (!_stack.empty()) {
		/* get next item from stack */
		_Item & item = *(_stack.top());

		/* valiate constraint */
		if (item.sc_state.size() < limit) {
			pop_stat();
		}
		/* reach the end of the children */
		else if (item.cursor >= item.op_state.size()) {
			/* find better solution */
			if (item.value >= maxvalue) {
				maxvalue = item.value;
				gen_re_solution(SufSOP);
			}
			pop_stat();	/* roll back */
		}
		/* otherwise, search further */
		else {
			rpsh_stat(); item.cursor = item.cursor + 1;
		}

	} /* end while : stack solution */

	  /* validate solution */
	if (SufSOP.empty()) {
		CError error(CErrorType::Runtime,
			"SuOperatorSearch::suf_solve",
			"Solution is not found");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else return;
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
void SuOperatorSearch::init_stat(const std::set<std::string> & init) {
	/* clear and initialize state */
	while (!_stack.empty()) pop_stat();
	_Item * item = new _Item();

	/* for operator-candidates */
	auto obeg = init.begin(), oend = init.end();
	while (obeg != oend) {
		const std::string & op = *(obeg++);
		if (op_scs.count(op) == 0)
			continue;
		else
			item->op_state.push_back(op);
	}
	/* compute covered clusters */
	for (int i = 0; i < item->op_state.size(); i++) {
		const std::string & op = item->op_state[i];
		if (op_scs.count(op) > 0) {
			auto op_iter = op_scs.find(op);
			const std::set<MuCluster *> & scs = *(op_iter->second);
			auto scs_beg = scs.begin(), scs_end = scs.end();
			while (scs_beg != scs_end) {
				MuCluster * covelm = *(scs_beg++);
				item->sc_state.insert(covelm);
			}
		}
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
void SuOperatorSearch::rpsh_stat() {
	if (_stack.empty()) {
		CError error(CErrorType::Runtime,
			"SuOperatorSearch::rpsh_stat",
			"Incorrect access: empty state");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else {
		/* get next operator to be selected */
		_Item & parent = *(_stack.top());

		/* create new item into stack */
		_Item & child = *(new _Item());

		/* compute operator and clusters */
		for (int i = 0; i < parent.op_state.size(); i++) {
			if (i != parent.cursor) {
				/* insert next operator from parent */
				const std::string & op = parent.op_state[i];
				child.op_state.push_back(op);

				/* update subsuming clusters */
				auto iter = op_scs.find(op);
				const std::set<MuCluster *> & clusters = *(iter->second);
				auto beg = clusters.begin(), end = clusters.end();
				while (beg != end) 
					child.sc_state.insert(*(beg++));
			}
		} /* end for op-sc */

		/* to ignore previous used operator */
		child.cursor = parent.cursor;

		/* compute value */
		const std::string & rop = parent.op_state[parent.cursor];
		child.value = parent.value + eval(rop);

		/* push in last */ _stack.push(&child);
	}
}
void SuOperatorSearch::gpsh_stat() {
	if (_stack.empty()) {
		CError error(CErrorType::Runtime,
			"SuOperatorSearch::rpsh_stat",
			"Incorrect access: empty state");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else {
		/* get next operator to be selected */
		_Item & parent = *(_stack.top());

		/* create the next item to be pushed */
		_Item & child = *(new _Item());

		/* get the selected operator */
		if (parent.cursor >= parent.op_state.size()) {
			CError error(CErrorType::OutOfIndex,
				"SuOperator::gpsh_stat()",
				"Out of limit for parent.op_stat");
			CErrorConsumer::consume(error);
			exit(CErrorType::OutOfIndex);
		}
		const std::string & op = parent.op_state[parent.cursor];

		/* compute subsuming clusters */
		auto op_iter = op_scs.find(op); child.sc_state.clear();
		const std::set<MuCluster *> & reduced = *(op_iter->second);
		auto sc_beg = parent.sc_state.begin();
		auto sc_end = parent.sc_state.end();
		while (sc_beg != sc_end) {
			MuCluster * cluster = *(sc_beg++);
			if (reduced.count(cluster) == 0)
				child.sc_state.insert(cluster);
		}

		/* compute operators */
		if (!child.sc_state.empty()) {
			auto next_iter = child.sc_state.begin();
			MuCluster * first = *(next_iter);
			auto first_iter = sc_ops.find(first); child.op_state.clear();
			const std::set<std::string> & first_ops = *(first_iter->second);
			auto first_beg = first_ops.begin(), first_end = first_ops.end();
			while (first_beg != first_end) 
				child.op_state.push_back(*(first_beg++));
		}

		/* compute cursor and value */
		child.value = parent.value + eval(op);
		child.cursor = 0;

		/* push child */ _stack.push(&child);
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
void SuOperatorSearch::gen_re_solution(std::set<std::string> & SOP) {
	SOP.clear(); _Item & top = *(_stack.top());
	for (int i = 0; i < top.op_state.size(); i++) 
		SOP.insert(top.op_state[i]);
}
void SuOperatorSearch::gen_coverset(const std::set<std::string> & ops, std::set<MuCluster *> & clusters) {
	clusters.clear();

	auto op_beg = ops.begin(), op_end = ops.end();
	while (op_beg != op_end) {
		const std::string & op = *(op_beg++);
		if (op_scs.count(op) > 0) {
			auto op_iter = op_scs.find(op);
			const std::set<MuCluster *> & scs = *(op_iter->second);
			auto scs_beg = scs.begin(), scs_end = scs.end();
			while (scs_beg != scs_end) {
				MuCluster * covelm = *(scs_beg++);
				clusters.insert(covelm);
			}
		}
	}
}
double SuOperatorSearch::get_coverage(const std::set<std::string> & ops) {
	std::set<MuCluster *> clusters;
	auto op_beg = ops.begin(), op_end = ops.end();
	while (op_beg != op_end) {
		const std::string & op = *(op_beg++);
		if (op_scs.count(op) > 0) {
			auto op_iter = op_scs.find(op);
			const std::set<MuCluster *> & scs = *(op_iter->second);
			auto scs_beg = scs.begin(), scs_end = scs.end();
			while (scs_beg != scs_end) {
				MuCluster * covelm = *(scs_beg++);
				clusters.insert(covelm);
			}
		}
	}

	size_t covernum = clusters.size();
	size_t totalnum = context.get_global_group().get_subsumings().size();
	return ((double)covernum) / ((double)totalnum);
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
	// eliminate by now (for print dominator score map)
//	std::ofstream out3(dir->get_path() + "/opscores.txt");
//	std::ofstream out4(dir->get_path() + "/avgscores.txt");
//	this->write_contr_domscore(data, SOP, out3, out4); 
//	out3.close(); out4.close();
//	std::ofstream out4(dir->get_path() + "/ctscores.txt");
//	this->write_domscore_line(data, out4); out4.close();
//	std::ofstream out5(dir->get_path() + "/mutants.txt");
//	this->write_mutants(data, out5); out5.close();

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
void OperatorWriter::write_contr_domscore(MutLevel & context, 
	const std::set<std::string> & SOP, std::ostream & out, std::ostream & aout) {

	/* initialization */
	SuOperatorSearch search(context); search.start();
	TestMachine tmachine(context); BitSeq seq(SOP.size());
	std::set<std::string> sel_ops; 
	TestSet * testset = ctest.malloc_test_set();

	std::set<MuCluster *> coverset;
	std::map<unsigned, std::vector<double> *> solutions;

	/* compute all pairs for solutions */
	unsigned long long limit = std::pow(2, SOP.size());
	for (unsigned long long k = 0; k < limit; k++) {
		/* generate selected operators */
		gen_combination(SOP, seq, sel_ops);	

		/* get cover-set and numbers */
		search.gen_coverset(sel_ops, coverset);	
		unsigned covernum = coverset.size();

		/* compute dominator score */
		tmachine.generate_by_operators(*testset, sel_ops);
		double score = tmachine.evaluate(*testset);
		
		/* get the set of scores for the coverage */
		if (solutions.count(covernum) == 0) 
			solutions[covernum] = new std::vector<double>();
		auto iter = solutions.find(covernum);
		std::vector<double> & scores = *(iter->second);

		/* add score */ scores.push_back(score);
		/* roll to the next */ seq.increase();	

	} /* end for */

	/* clear resources */
	search.finish(); ctest.delete_test_set(testset);

	/* title */ 
	// out << "contribution\tmin\tavg\tmax\n";
	out << "contribution\tdominator score\n";

	/* output solution */
	unsigned int totalnum = context.get_global_group().get_subsumings().size();
	auto sbeg = solutions.begin(), send = solutions.end();
	std::map<double, double> avg_solutions;
	while (sbeg != send) {
		/* get scores and coverage */
		unsigned covernum = sbeg->first;
		std::vector<double> * scores = sbeg->second;

		/* get coverage */
		double coverage = ((double)covernum) / ((double)totalnum);

		/* get score set for output */
		auto scbeg = scores->begin(), scend = scores->end();
		double average = 0;
		while (scbeg != scend) {
			double score = *(scbeg++);
			out << coverage << "\t" << score << "\n";
			average += score;
		}
		average = average / scores->size();
		avg_solutions[coverage] = average;

		/*
		double min = 1.0, max = 0.0, sum = 0.0, avg;
		auto scbeg = scores->begin(), scend = scores->end();
		while (scbeg != scend) {
			double score = *(scbeg++);
			if (min > score) min = score;
			if (max < score) max = score;
			sum += score;
		}
		avg = sum / scores->size();

		out << coverage << "\t" << min << "\t" << avg << "\t" << max << "\n";
		*/

		delete scores;	/* delete resource */
		sbeg++;	/* roll to the next */
	}
	solutions.clear();

	/* output average */
	aout << "contribution\taverage score\n";
	auto abeg = avg_solutions.begin();
	auto aend = avg_solutions.end();
	while (abeg != aend) {
		double coverage = abeg->first;
		double score = abeg->second;
		aout << coverage << "\t" << score << "\n";
		abeg++;
	}

	/* end */ return;
}
void OperatorWriter::gen_combination(const std::set<std::string> & source,
	const BitSeq & sel_bits, std::set<std::string> & target) {

	/* initialization */
	target.clear(); BitSeq::size_t k = 0;
	auto beg = source.begin(), end = source.end();

	/* select operators */
	while (beg != end) {
		const std::string & op = *(beg++);
		if (sel_bits.get_bit(k++) == BIT_1) 
			target.insert(op);
	}
}
void OperatorWriter::write_domscore_line(MutLevel & context, std::ostream & out) {
	/* initialization */
	const std::set<MuCluster *> requirements 
		= context.get_global_group().get_subsumings();
	std::vector<MuCluster *> requirement_list;
	auto rbeg = requirements.begin(), rend = requirements.end();
	while (rbeg != rend) requirement_list.push_back(*(rbeg++));

	/* solutions for contribution to dominator score(s) */
	std::map<BitSeq::size_t, std::set<double> *> solutions;
	std::set<MuCluster *> sel_reqs; TestMachine tmachine(context);
	TestSet * tests = ctest.malloc_test_set();

	/* combination iterate */
	BitSeq seq(requirements.size());
	do {
		/* generate next requirement for selection */
		gen_combination(requirement_list, seq, sel_reqs);	

		/* evaluate test dominator score */
		tmachine.generate_by_requirement(*tests, sel_reqs);
		double dom_score = tmachine.evaluate(*tests);

		/* get the set of dominator scores */
		BitSeq::size_t req_num = seq.degree();
		if (solutions.count(req_num) == 0) 
			solutions[req_num] = new std::set<double>();
		auto iter = solutions.find(req_num);
		std::set<double> & scores = *(iter->second);

		/* insert score */ scores.insert(dom_score);

		/* roll to the next */ seq.increase();
	} while (!seq.all_zeros());

	/* delete test set */ ctest.delete_test_set(tests);

	/* title */ out << "contribution\tmin\tavg\tmax\n";

	/* output solution */
	unsigned int totalnum = requirements.size();
	auto sbeg = solutions.begin();
	auto send = solutions.end();
	while (sbeg != send) {
		/* get scores and coverage */
		unsigned covernum = sbeg->first;
		std::set<double> * scores = sbeg->second;
		double coverage = ((double)covernum) / ((double)totalnum);

		double min = 1.0, max = 0.0, sum = 0.0, avg;
		auto scbeg = scores->begin(), scend = scores->end();
		while (scbeg != scend) {
			double score = *(scbeg++);
			if (min > score) min = score;
			if (max < score) max = score;
			sum += score;
		}
		avg = sum / scores->size();

		out << coverage << "\t" << min << "\t" << avg << "\t" << max << "\n";

		delete scores;	/* delete resource */
		sbeg++;	/* roll to the next */
	}
	solutions.clear();
}
void OperatorWriter::gen_combination(const std::vector<MuCluster *> & source,
	const BitSeq & sel_bits, std::set<MuCluster *> & target) {
	/* initialization */
	target.clear(); BitSeq::size_t k = 0;
	auto beg = source.begin(), end = source.end();

	/* select operators */
	while (beg != end) {
		MuCluster * cluster = *(beg++);
		if (sel_bits.get_bit(k++) == BIT_1)
			target.insert(cluster);
	}
}
void OperatorWriter::trim_spaces(std::string & text) {
	std::string cache; char ch; int n = text.length();
	for (int i = 0; i < n; i++) {
		ch = text[i];

		if (ch == '\t' || ch == '\n') continue;
		else cache += ch;
	}
	text = cache;
}
void OperatorWriter::write_mutants(MutLevel & data, std::ostream & out) {
	/* declaration */
	MutantSpace & mspace = data.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	const MutGroup & group = data.get_global_group();

	out << "id\toperator\torigin\treplace\tcategory\tmode\n";
	for (mid = 0; mid < num; mid++) {
		/* get next mutant */
		Mutant & mutant = mspace.get_mutant(mid);
		const std::string & op = mutant.get_operator();
		if (!data.has_operator(op)) continue;
		const MutGroup & op_group = data.get_operator_group(op);

		/* get category */
		std::string category;
		if (group.category_of(mid) == Subsuming_Category) {
			category = "Subsuming";
		}
		else if (op_group.category_of(mid) == Equvalent_Category) {
			category = "Equivalent";
		}
		else continue;

		/* text */
		const Mutation & mutation = 
			mutant.get_mutation(mutant.get_orders() - 1);
		std::string origin = mutation.get_location().get_text_at();
		std::string replace = mutation.get_replacement();
		trim_spaces(origin); trim_spaces(replace);

		/* output line */
		out << mid << "\t" << op << "\t";
		out << origin << "\t" << replace << "\t";
		out << category << "\t" << "";
		out << "\n";
	}

	/* return */ out << std::endl;
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
	//std::cerr << "Dominator score: " << K << "/" << M << "\n";

	/* return */ return ((double) K) / ((double) M);
}
bool TestMachine::is_killed(const TestSet & tests, const MuCluster & cluster) {
	const BitSeq & tseq = tests.get_set_vector();
	const BitSeq & mseq = cluster.get_score_vector();
	BitSeq rseq(mseq); rseq.conjunct(tseq);
	//std::cerr << "\tResult: " << rseq.degree() << "\t" << rseq.all_zeros() << "\t" << !(rseq.all_zeros()) << "\n";
	return !(rseq.all_zeros());
}
void TestMachine::generate_by_operators(TestSet & tests, const std::set<std::string> & ops) {
	/* selected operators */
	std::set<std::string> operators;
	auto beg = ops.begin(), end = ops.end();
	while (beg != end) {
		const std::string & op = *(beg++);
		if (context.has_operator(op))
			operators.insert(op);
	}

	/* add all subsuming clusters into requirements */
	std::set<MuCluster *> requirements;
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

	/* generate tests by requirements */
	generate_by_requirement(tests, requirements);

	/* return */ return;
}
void TestMachine::generate_by_requirement(TestSet & tests, const std::set<MuCluster *> & requirements) {
	this->greedy_generate_tests(tests, requirements);
}
void TestMachine::greedy_generate_tests(TestSet & tests, const std::set<MuCluster *> & RS) {
	/* initialization */
	tests.clear(); std::set<MuCluster *> reduced;
	std::set<MuCluster *> requirements;
	auto rbeg = RS.begin(), rend = RS.end();
	while (rbeg != rend) 
		requirements.insert(*(rbeg++));

	while (!requirements.empty()) {
		/* get the next requirement to be met */
		MuCluster & next_req = *(*(requirements.begin()));

		/* eliminate equivalent cluster */
		if (next_req.get_score_degree() == 0) {
			requirements.erase(&next_req);
			continue;
		}

		/* get the test set that kill this requirement */
		const BitSeq & score_set = next_req.get_score_vector();
		BitSeq::size_t score_deg = next_req.get_score_degree();

		/* get one test that kill this requirement */
		BitSeq::size_t rand_seed = gen_random_seed(score_deg);
		TestCase::ID tid = find_test_at(score_set, rand_seed);

		/* update requirement with newly test */
		reduced.clear();
		auto beg = requirements.begin();
		auto end = requirements.end();
		while (beg != end) {
			MuCluster * req = *(beg++);
			const BitSeq & tset = req->get_score_vector();
			if (tset.get_bit(tid) == BIT_1) 
				reduced.insert(req);
		}

		/* eliminate killed requirements */
		beg = reduced.begin(), end = reduced.end();
		while (beg != end) 
			requirements.erase(*(beg++));

		/* add into test set */ tests.add_test(tid);

	} /* end while: kill requirements */

	/* return */ return;
}
BitSeq::size_t TestMachine::gen_random_seed(BitSeq::size_t n) {
	return 0 + std::rand() % n;
}
TestCase::ID TestMachine::find_test_at(const BitSeq & bits, TestCase::ID seed) {
	TestCase::ID k, n = bits.bit_number(); 
	unsigned origin_seed = seed;
	for (k = 0; k < n; k++) {
		if (bits.get_bit(k) == BIT_1) {
			if (seed == 0) return k;
			else seed--;
		}
	}

	/* not found */
	CError error(CErrorType::InvalidArguments, "TestMachine::find_test_at", 
		"Undefined index (" + std::to_string(origin_seed) + ")");
	CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
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
/* identify and output the alpha-suf SMOs 
		{alpha, contribution, dominator score, mutants, operators}
*/
static void output_operators(
	const MutLevel & data, SuOperatorSearch & op_search,
	const std::set<std::string> & ops, std::ostream & out,
	CTest & ctest) {
	/* get the coverage of killed subsuming clusters */
	double coverage = op_search.get_coverage(ops);

	/* compute dominator score */
	TestMachine tmachine(data);
	TestSet * tests = ctest.malloc_test_set();
	tmachine.generate_by_operators(*tests, ops);
	double score = tmachine.evaluate(*tests);
	ctest.delete_test_set(tests);

	/* get the number of mutants */
	auto beg = ops.begin(), end = ops.end();
	out << "operators: \t{"; size_t num = 0;
	while (beg != end) {
		const std::string & op = *(beg++);
		if (data.has_operator(op)) {
			out << op << ";";
			const MutGroup & group = data.get_operator_group(op);
			num += group.get_mutants().number_of_mutants();
		}
	}
	out << "}\n";
	size_t total = data.get_space().number_of_mutants();

	/* output */
	out << "mutants: \t" << ((double) num) / ((double) total)  << "\n";
	out << "contribution: \t" << coverage << "\n";
	out << "dominator score; \t" << score << "\n";

	/* return */ return;
}
/* get E-selective operators */
static void E_selective_operators(std::set<std::string> & eops) {
	eops.insert("u-OAAN");
	eops.insert("u-ORRN");
	eops.insert("u-OLLN");
	eops.insert("I-DirVarIncDec");
	// eops.insert("I-IndVarIncDec");
}

/* test */
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "replace";
	TestType ttype = TestType::replace;

	// get root file and analysis dir 
	File & root = *(new File(prefix + prname));
	 
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

		// debug operators
		auto cbeg = debug_operators.begin();
		auto cend = debug_operators.end();
		std::cout << "\nUnused operators (" << debug_operators.size() << ")\n";
		while (cbeg != cend) {
			const std::string & op = *(cbeg++);
			std::cout << op << "\n";
		}
		std::cout << "\n";

		// solve subsuming operators
		std::set<std::string> SOP;
		SuOperatorSearch search(data);
		search.start(); 
		// search.solve(SOP); 
		// search.re_solve(SOP);
		search.gr_solve(SOP);
		search.finish();
		std::cout << "Solving subsuming operators: " << SOP.size() << "\n";

		// end output
		OperatorWriter writer(ctest);
		writer.open(root);
		writer.write(data, SOP);
		writer.close();
		std::cout << "Output subsuming mutants: finished\n";

		// output sufficient operators, SOP, and EOPs
		std::ofstream out(root.get_path() + "/analysis/sufsop.txt");

		search.start();
		out << "alpha: 100%\n";
		output_operators(data, search, SOP, out, ctest);
		std::cout << "alpha 100% finished (" << SOP.size() << ")\n";
		out << "\n";

		std::set<std::string> SufSOP;

		search.suf_solve(SOP, SufSOP, 0.60);
		std::cout << "Compute 0.6-sufficient SOPs...\n";

		out << "alpha: 60%\n";
		output_operators(data, search, SufSOP, out, ctest);
		out << "\n";

		search.suf_solve(SOP, SufSOP, 0.80);
		std::cout << "Compute 0.8-sufficient SOPs...\n";

		out << "alpha: 90%\n";
		output_operators(data, search, SufSOP, out, ctest);
		out << "\n";

		search.suf_solve(SOP, SufSOP, 0.90);
		std::cout << "Compute 0.9-sufficient SOPs...\n";

		out << "alpha: 90%\n";
		output_operators(data, search, SufSOP, out, ctest);
		out << "\n";

		std::set<std::string> eops;
		E_selective_operators(eops);
		out << "E-selective\n";
		output_operators(data, search, eops, out, ctest);
		out << "\n";

		search.finish(); out.close();

		std::cout << "Output finished.\n\n";
	}

	// delete memory
	delete &cscore;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	// exit
	std::cout << "Press any key to exit...\n"; getchar(); exit(0);
}
