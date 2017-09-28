#include "suoprt.h"

MuClusterSet::MuClusterSet(const MSGraph & g) : graph(g) {
	/* compute equivalent cluster */
	eq_cluster = nullptr;
	if (g.get_roots().size() == 1) {
		auto beg = g.get_roots().begin();
		MuCluster * cluster = *(beg);
		if (cluster->get_score_vector().all_zeros())
			eq_cluster = cluster;
	}

	/* compute subsuming clusters */
	subsuming_clusters.clear();
	if (eq_cluster == nullptr) {
		auto beg = g.get_roots().begin();
		auto end = g.get_roots().end();
		while (beg != end) 
			subsuming_clusters.insert(*(beg++));
	}
	else {
		const std::vector<MuSubsume> & edges = eq_cluster->get_ou_port().get_edges();
		for (int i = 0; i < edges.size(); i++) {
			const MuSubsume & edge = edges[i];
			MuCluster * target = &(edge.get_target());
			subsuming_clusters.insert(target);
		}
	}

	/* compute subsumed clusters */
	subsumed_clusters.clear();
	int len = g.size();
	for (int i = 0; i < len; i++) {
		MuCluster * cluster = &(g.get_cluster(i));
		if (cluster != eq_cluster &&
			subsuming_clusters.count(cluster) == 0) {
			subsumed_clusters.insert(cluster);
		}
	}
}
MuCluster & MuClusterSet::get_equivalents() const {
	if (eq_cluster == nullptr) {
		CError error(CErrorType::InvalidArguments, 
			"MuClusterSet::get_equivalents", 
			"Invalid access: no equivalents");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else return *eq_cluster;
}
char MuClusterSet::category_of(MuCluster & cluster) const {
	MuCluster * cptr = &cluster;
	if (cptr == eq_cluster) return Equivalent;
	else if (subsuming_clusters.count(cptr) > 0) return Subsuming;
	else if (subsumed_clusters.count(cptr) > 0) return Subsumed;
	else return NotBelong;
}

OpClusterMap::OpClusterMap(const MSGraph & g) : graph(g) {
	MutantSpace & mspace = graph.get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();

	for (mid = 0; mid < num; mid++) {
		Mutant & mutant = mspace.get_mutant(mid);			/* get next mutant */
		const std::string & op = mutant.get_operator();		/* get its operator */

		if (operators.count(op) == 0) {
			operators.insert(op);
			opc_map[op] = new std::set<MuCluster *>();
		}

		if (graph.has_cluster_of(mid)) {
			MuCluster * cluster = &(graph.get_cluster_of(mid));

			if (clusters.count(cluster) == 0) {
				clusters.insert(cluster);
				cop_map[cluster] = new std::set<std::string>();
			}

			auto op_iter = opc_map.find(op);
			std::set<MuCluster *> & cset = *(op_iter->second);
			cset.insert(cluster);

			auto cl_iter = cop_map.find(cluster);
			std::set<std::string> & oset = *(cl_iter->second);
			oset.insert(op);
		}
	}
}
OpClusterMap::~OpClusterMap() {
	auto beg1 = cop_map.begin();
	auto end1 = cop_map.end();
	while (beg1 != end1) 
		delete (beg1++)->second;

	auto beg2 = opc_map.begin();
	auto end2 = opc_map.end();
	while (beg2 != end2) 
		delete (beg2++)->second;

	clusters.clear();
	operators.clear();
	cop_map.clear();
	opc_map.clear();
}
const std::set<MuCluster *> & OpClusterMap::get_clusters_of(const std::string & op) const {
	if (opc_map.count(op) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"OpClusterMap::get_clusters_of", 
			"Invalid operator: \"" + op + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = opc_map.find(op);
		return *(iter->second);
	}
}
const std::set<std::string> & OpClusterMap::get_operators_of(MuCluster & cluster) const {
	if (cop_map.count(&cluster) == 0) {
		CError error(CErrorType::InvalidArguments,
			"OpClusterMap::get_operators_of",
			"Invalid clusters");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = cop_map.find(&cluster);
		return *(iter->second);
	}
}

OpMutantMap::OpMutantMap(MutantSpace & space) : mspace(space) {
	Mutant::ID mid, num = mspace.number_of_mutants();
	for (mid = 0; mid < num; mid++) {
		Mutant & mutant = space.get_mutant(mid);
		const std::string & op = mutant.get_operator();

		if (operators.count(op) == 0) {
			operators.insert(op);
			op_map[op] = new std::set<Mutant::ID>();
		}

		auto iter = op_map.find(op);
		std::set<Mutant::ID> & ids = *(iter->second);
		ids.insert(mid);
	}
}
const std::set<Mutant::ID> & OpMutantMap::get_mutants_of(const std::string & op) const {
	if (op_map.count(op) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"OpMutantMap::get_mutants_of", 
			"Invalid operator: \"" + op + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = op_map.find(op);
		return *(iter->second);
	}
}

SOperatorSet::SOperatorSet(
	const MuClusterSet & cset, const OpClusterMap & map, const OpMutantMap & ops) 
	: clusters(cset), mappings(map), operators(ops) {
	/* compute SOPSet */
	solve_by_greedy(1.0, SOPSet);
}
SOperatorSet::~SOperatorSet() {
	SOPSet.clear();
	while (!_stack.empty()) {
		delete _stack.top();
		_stack.pop();
	}
}
bool SOperatorSet::solve_by_bounds(
	double alpha, std::set<std::string> & target) {
	/* compute the constraint */
	size_t limit = clusters.get_subsuming().size();
	size_t bound = limit * alpha + 1;
	if (bound > limit) bound = limit;
	unsigned reduced_num = 0; target.clear();

	/* get initial seed */
	std::set<std::string> seed;
	auto beg = clusters.get_subsuming().begin();
	auto end = clusters.get_subsuming().end();
	while (beg != end) {
		MuCluster * cluster = *(beg++);
		const std::set<std::string> & ops =
			mappings.get_operators_of(*cluster);
		auto obeg = ops.begin(), oend = ops.end();
		while (obeg != oend)
			seed.insert(*(obeg++));
	}

	/* compute solutions */
	bound_init_item(seed);	/* initialize stack */
	while (!_stack.empty()) {
		_Item & top = *(_stack.top());	/* get the parent item */

		/* invalid constraints */
		if (top.sc_state.size() < bound) {
			pope_item();
		}
		/* finish all its children */
		else if (top.cursor >= top.op_state.size()) {
			if (top.value >= reduced_num) {
				reduced_num = top.value;
				top_solutions(target);
			}
			pope_item();
		}
		/* find deep-search */
		else {
			bound_push_item(); 
			top.cursor = top.cursor + 1;
		}

	} /* end while */

	/* validate solution */
	if (target.empty()) {
		CError error(CErrorType::Runtime, 
			"SOperatorSet::solve_by_bounds", 
			"Internal errors");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else return true;
}
bool SOperatorSet::solve_by_greedy(
	double alpha, std::set<std::string> & target) {
	/* compute the constraint */
	unsigned bound = mappings.get_graph().
		get_space().number_of_mutants();

	/* get the lower-bound of constraints */
	int limit = clusters.get_subsuming().size();
	limit = limit * (1 - alpha);
	if (limit < 0) limit = 0;

	/* initialize solutions */
	target.clear(); 

	/* compute solutions */
	greedy_init_item(clusters.get_subsuming());
	while (!_stack.empty()) {
		_Item & parent = *(_stack.top());	/* get parent item */

		/* validate bounds: for number of mutants */
		if (parent.value > bound) {
			pope_item();
		}
		/* find feasible solution */
		else if (parent.sc_state.size() <= limit) {
			bound = parent.value;
			pope_item();
			lin_solutions(target);
		}
		/* finish all its children */
		else if (parent.cursor >= parent.op_state.size()) {
			pope_item();
		}
		/* further deep search */
		else {
			greedy_push_item();
			parent.cursor = parent.cursor + 1;
		}
	}

	/* validate solution */
	if (target.empty()) {
		CError error(CErrorType::Runtime,
			"SOperatorSet::solve_by_greedy",
			"Internal errors");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else return true; 
}
unsigned SOperatorSet::evaluate(const std::string & op) const {
	if (operators.has_operator(op)) 
		return operators.get_mutants_of(op).size();
	else return 0;
}
void SOperatorSet::pope_item() {
	if (_stack.empty()) {
		CError error(CErrorType::Runtime, 
			"SOperatorSet::pope_item", 
			"_stack is empty");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else {
		delete _stack.top();
		_stack.pop();
	}
}
void SOperatorSet::bound_init_item(const std::set<std::string> & seed) {
	/* clear stack */
	while (!_stack.empty()) {
		delete _stack.top();
		_stack.pop();
	}

	/* create a new root */
	_Item & root = *(new _Item());	

	/* update operators */
	auto sbeg = seed.begin();
	auto send = seed.end();
	while (sbeg != send)
		root.op_state.push_back(*(sbeg++));

	/* compute initial requirement */
	std::set<MuCluster *> space = clusters.get_subsuming();
	for (int i = 0; i < root.op_state.size(); i++) {
		/* get the next operator and insert to child solution space */
		const std::string & opi = root.op_state[i];

		/* compute covered subsuming clusters */
		if (!space.empty()) {
			/* get the subsuming clusters covered by this operator */
			auto beg = space.begin(), end = space.end();
			while (beg != end) {
				MuCluster * cluster = *(beg++);
				const std::set<std::string> & ops
					= mappings.get_operators_of(*cluster);
				if (ops.count(opi) > 0)
					root.sc_state.insert(cluster);
			}

			/* eliminate covered clusters */
			beg = root.sc_state.begin();
			end = root.sc_state.end();
			while (beg != end)
				space.erase(*(beg++));
		}
	}	/* end for */

	/* compute value */
	root.value = 0; root.cursor = 0;

	/* push to the stack */
	_stack.push(&root);	
}
void SOperatorSet::bound_push_item() {
	/* validate 1: no parent */
	if (_stack.empty()) {
		CError error(CErrorType::Runtime, 
			"SOperatorSet::bound_push_item", 
			"stack is empty");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	/* validate 2: operator selected */
	_Item & parent = *(_stack.top());
	if (parent.cursor >= parent.op_state.size()) {
		CError error(CErrorType::OutOfIndex,
			"SOperatorSet::bound_push_item",
			"parent's cursor is out of index");
		CErrorConsumer::consume(error);
		exit(CErrorType::OutOfIndex);
	}

	/* get selected operator (for elimination) */
	const std::string & op = parent.op_state[parent.cursor];
	std::set<MuCluster *> space = clusters.get_subsuming();

	/* get new child */ 
	_Item & child = *(new _Item()); child.cursor = 0;

	/* compute value */
	child.value = parent.value + evaluate(op);

	/* compute operator set and subsuming clusters */
	for (int i = parent.cursor + 1; i < parent.op_state.size(); i++) {
		/* get the next operator and insert to child solution space */
		const std::string & opi = parent.op_state[i];
		child.op_state.push_back(opi);

		/* compute covered subsuming clusters */
		if (!space.empty()) {
			/* get the subsuming clusters covered by this operator */
			auto beg = space.begin(), end = space.end();
			while (beg != end) {
				MuCluster * cluster = *(beg++);
				const std::set<std::string> & ops
					= mappings.get_operators_of(*cluster);
				if (ops.count(opi) > 0)
					child.sc_state.insert(cluster);
			}

			/* eliminate covered clusters */
			beg = child.sc_state.begin();
			end = child.sc_state.end();
			while (beg != end)
				space.erase(*(beg++));
		}
	} /* end for */

	/* push the stack */ _stack.push(&child);
}
void SOperatorSet::top_solutions(std::set<std::string> & target) {
	if (_stack.empty()) {
		CError error(CErrorType::InvalidArguments,
			"SOperatorSet::top_solutions",
			"Invalid access: stack is empty");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		_Item & top = *(_stack.top()); target.clear();
		for (int i = 0; i < top.op_state.size(); i++)
			target.insert(top.op_state[i]);
	}
}
void SOperatorSet::greedy_init_item(const std::set<MuCluster *> & seed) {
	/* clear stack */
	while (!_stack.empty()) {
		delete _stack.top();
		_stack.pop();
	}

	/* create root */
	_Item & root = *(new _Item()); root.cursor = 0;

	/* compute clusters */
	auto beg = seed.begin(), end = seed.end();
	while (beg != end) {
		MuCluster * cluster = *(beg++);
		if (clusters.category_of(*cluster) 
			== MuClusterSet::Subsuming) 
			root.sc_state.insert(cluster);
	}

	/* compute operators for selection */
	if (!root.sc_state.empty()) {
		MuCluster * first = *(root.sc_state.begin());
		const std::set<std::string> & ops 
			= mappings.get_operators_of(*first);
		auto obeg = ops.begin(), oend = ops.end();
		while (obeg != oend) 
			root.op_state.push_back(*(obeg++));
	}

	/* compute values (none is selected yet) */
	root.value = 0;

	/* push to stack */ _stack.push(&root);
}
void SOperatorSet::greedy_push_item() {
	/* validate 1: no parent */
	if (_stack.empty()) {
		CError error(CErrorType::Runtime,
			"SOperatorSet::bound_push_item",
			"stack is empty");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	/* validate 2: operator selected */
	_Item & parent = *(_stack.top());
	if (parent.cursor >= parent.op_state.size()) {
		CError error(CErrorType::OutOfIndex,
			"SOperatorSet::bound_push_item",
			"parent's cursor is out of index");
		CErrorConsumer::consume(error);
		exit(CErrorType::OutOfIndex);
	}

	/* get selected operator */
	const std::string & op = parent.op_state[parent.cursor];
	const std::set<MuCluster *> & covers = mappings.get_clusters_of(op);

	/* create new item */
	_Item & child = *(new _Item()); child.cursor = 0;

	/* update child.sc_state */
	auto beg = parent.sc_state.begin();
	auto end = parent.sc_state.end();
	while (beg != end) {
		MuCluster * cluster = *(beg++);
		if (covers.count(cluster) == 0)
			child.sc_state.insert(cluster);
	}

	/* update child.op_state */
	if (!child.sc_state.empty()) {
		MuCluster * first = *(child.sc_state.begin());
		const std::set<std::string> & ops 
			= mappings.get_operators_of(*first);
		auto obeg = ops.begin(), oend = ops.end();
		while (obeg != oend)
			child.op_state.push_back(*(obeg++));
	}

	/* get the number of selected mutants */
	child.value = parent.value + evaluate(op);

	/* push to stack */ _stack.push(&child);
}
void SOperatorSet::lin_solutions(std::set<std::string> & target) {
	target.clear();	/* clear the solution */

	std::stack<_Item *> cache;
	while (!_stack.empty()) {
		/* get the item of stack and cache */
		_Item * item = _stack.top();
		_stack.pop(); cache.push(item);

		/* get selected operator */
		const std::string & op = 
			item->op_state[item->cursor - 1];
		target.insert(op);
	}

	/* recover stack */
	while (!cache.empty()) {
		_Item * item = cache.top();
		cache.pop(); _stack.push(item);
	}
}

void SuOprtWriter::write(SOperatorSet & data) {
	if (dir == nullptr) {
		CError error(CErrorType::Runtime, 
			"SuOprtWriter::write", 
			"File is not opened");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}
	else {
		std::ofstream out1(dir->get_path() + "/summary.txt");
		write_summary(data, out1); out1.close();
		std::ofstream out2(dir->get_path() + "/coverage.txt");
		write_coverage(data, out2); out2.close();
		//std::ofstream out3(dir->get_path() + "/score_lines.txt");
		//write_scoreln(data, out3); out3.close();
		std::ofstream out4(dir->get_path() + "/op_scores.txt");
		write_op_scores(data, out4); out4.close();
		std::ofstream out5(dir->get_path() + "/mutantset.txt");
		write_mutants(data, out5); out5.close();
	}
}
void SuOprtWriter::write_summary(SOperatorSet & data, std::ostream & out) {
	out << "#Mutants\t#Tests\t#Operators\t#SuMutants\t#SuOperators\t#SOpMutants\n";

	out << data.get_operators().get_space().number_of_mutants() << "\t";
	out << tspace.number_of_tests() << "\t";
	out << data.get_operators().get_operators().size() << "\t";
	out << data.get_clusters().get_subsuming().size() << "\t";
	out << data.get_subsuming_operators().size() << "\t";

	unsigned count = 0;
	auto beg = data.get_subsuming_operators().begin();
	auto end = data.get_subsuming_operators().end();
	while (beg != end) {
		const std::string & op = *(beg++);
		count += data.get_operators().get_mutants_of(op).size();
	}
	out << count << "\n";
}
void SuOprtWriter::write_coverage(SOperatorSet & data, std::ostream & out) {
	out << "operator\t#Mut\t#Equiv\t#Subsuming\t#Coverage\n";

	std::map<std::string, unsigned> equiv_map;
	MutantSpace & mspace = data.get_operators().get_space();
	Mutant::ID mid = 0, num = mspace.number_of_mutants();
	const MuClusterSet clusters = data.get_clusters();
	for (mid = 0; mid < num; mid++) {
		if (clusters.get_graph().has_cluster_of(mid)) {
			MuCluster & cluster = clusters.get_graph().get_cluster_of(mid);
			if (clusters.category_of(cluster) == MuClusterSet::Equivalent) {
				Mutant & mutant = mspace.get_mutant(mid);
				const std::string & op = mutant.get_operator();

				if (equiv_map.count(op) == 0)
					equiv_map[op] = 0;

				auto iter = equiv_map.find(op);
				unsigned size = iter->second;
				equiv_map[op] = size + 1;
			}
		}
	}

	auto beg = data.get_operators().get_operators().begin();
	auto end = data.get_operators().get_operators().end();
	while (beg != end) {
		const std::string & op = *(beg++);
		out << op << "\t";
		out << data.get_operators().get_mutants_of(op).size() << "\t";

		unsigned equivs = 0;
		if (equiv_map.count(op) > 0) {
			auto iter = equiv_map.find(op);
			equivs = iter->second;
		}
		out << equivs << "\t";

		std::set<MuCluster *> covered;
		get_coverage(data, covered, op);
		out << covered.size() << "\t";

		out << ((double)covered.size()) / ((double)clusters.get_subsuming().size()) << "\n";
	}	/* end while */
}
void SuOprtWriter::eval_operators(SOperatorSet & data, const std::vector<double> & alist) {
	/* validation */
	if (dir == nullptr) {
		CError error(CErrorType::Runtime, 
			"SuOprtWriter::eval_operators",
			"Invalid access: not opened");
		CErrorConsumer::consume(error);
		exit(CErrorType::Runtime);
	}

	/* declarations */
	std::set<std::string> suffOps; double alpha;
	TestSet * tests = ((CTest &)tspace.get_project()).malloc_test_set();
	std::ofstream out(dir->get_path() + "/operators.txt");

	/* title */
	out << "name\toperators\tmutants\tcoverage\tscore\n";

	for (int i = 0; i < alist.size(); i++) {
		double alpha = alist[i]; data.get_sufficient_operators_fast(alpha, suffOps);
		out << alpha << "-SMOs\t"; write_oplist(data, suffOps, out, *tests);
	}

	/* E-selective set */
	suffOps.clear(); 
	suffOps.insert("u-OAAN");
	suffOps.insert("u-OLLN");
	suffOps.insert("u-OLRN");
	suffOps.insert("u-ORRN");
	//suffOps.insert("I-DirVarIncDec");
	out << "E-selective\t";
	write_oplist(data, suffOps, out, *tests);


	/* close and return */ 
	((CTest &)tspace.get_project()).delete_test_set(tests);
	out.close();
}
void SuOprtWriter::write_oplist(SOperatorSet & data, 
	const std::set<std::string> & suffOps, 
	std::ostream & out, TestSet & tests) {
	std::set<MuCluster *> covers, cache;
	TestMachine tmachine(data); unsigned mutants = 0;
	auto beg = suffOps.begin(), end = suffOps.end();
	while (beg != end) {
		const std::string & op = *(beg++);
		if (!(data.get_operators().has_operator(op)))
			continue;

		out << op << "; ";

		get_coverage(data, cache, op);
		auto cbeg = cache.begin(), cend = cache.end();
		while (cbeg != cend) covers.insert(*(cbeg++));

		mutants += data.get_operators().get_mutants_of(op).size();
	}
	out << "\t";

	out << ((double) mutants) / ((double) data.get_operators().get_space().number_of_mutants()) << "\t";
	out << ((double)covers.size()) / ((double)data.get_clusters().get_subsuming().size()) << "\t";
	tmachine.generate_by_operators(tests, suffOps);
	out << tmachine.evaluate(tests) << "\n";
}
void SuOprtWriter::get_coverage(SOperatorSet & data, std::set<MuCluster *> & covers, const std::string & op) {
	covers.clear();

	const std::set<MuCluster *> & cset = 
		data.get_mappings().get_clusters_of(op);
	auto cbeg = cset.begin(), cend = cset.end();
	while (cbeg != cend) {
		MuCluster * cluster = *(cbeg++);
		if (data.get_clusters().category_of(
			*cluster) == MuClusterSet::Subsuming)
			covers.insert(cluster);
	}
}
void SuOprtWriter::write_mutants(SOperatorSet & data, std::ostream & out) {
	/* declarations */
	MutantSpace & mspace = data.get_operators().get_space();
	Mutant::ID mid, num = mspace.number_of_mutants();
	const MSGraph & graph = data.get_clusters().get_graph();

	/* title */
	out << "id\tcluster\toperator\tline\torigin\treplace\n";

	/* output lines */
	for (mid = 0; mid < num; mid++) {
		if (graph.has_cluster_of(mid)) {
			MuCluster & cluster = graph.get_cluster_of(mid);
			if (data.get_clusters().category_of(cluster) 
				== MuClusterSet::Subsuming) {
				/* get mutant information */
				Mutant & mutant = mspace.get_mutant(mid);
				const std::string & op = mutant.get_operator();
				const Mutation & mutation = mutant.
					get_mutation(mutant.get_orders() - 1);
				const CodeLocation & loc = mutation.get_location();

				std::string origin = loc.get_text_at();
				std::string replace = mutation.get_replacement();
				trim_spaces(origin); trim_spaces(replace);

				int line = loc.get_file().get_text()->lineOfIndex(loc.get_bias());

				out << mid << "\t" << cluster.get_id() << "\t" << op << "\t"
					<< line << "\t" << origin << "\t" << replace << "\n";
			}
		}
	}	/* end for */

	/* return */ return;
}
void SuOprtWriter::trim_spaces(std::string & text) {
	std::string cache; int n = text.length();
	for (int i = 0; i < n; i++) {
		char ch = text[i];
		if (ch == '\t' || ch == '\n')
			continue;
		else cache += ch;
	}
	text = cache;
}
void SuOprtWriter::gen_cov_score(SOperatorSet & data, 
	const std::set<std::string> & ops, std::ostream & out, TestSet & tests) {

	std::set<MuCluster *> covers, cache;
	TestMachine tmachine(data); 
	auto beg = ops.begin(), end = ops.end();
	while (beg != end) {
		const std::string & op = *(beg++);
		if (!(data.get_operators().has_operator(op)))
			continue;

		get_coverage(data, cache, op);
		auto cbeg = cache.begin(), cend = cache.end();
		while (cbeg != cend) covers.insert(*(cbeg++));
	}

	out << ((double)covers.size()) / ((double)data.get_clusters().get_subsuming().size()) << "\t";
	tmachine.generate_by_operators(tests, ops);
	out << tmachine.evaluate(tests) << "\n";
}
void SuOprtWriter::select_operators(const std::vector<std::string> & source, 
	const BitSeq & bits, std::set<std::string> & target) {
	target.clear();
	for (int i = 0; i < source.size(); i++) {
		if (bits.get_bit(i) == BIT_1) {
			target.insert(source[i]);
		}
	}
}
void SuOprtWriter::write_scoreln(SOperatorSet & data, std::ostream & out) {
	std::vector<std::string> oplist;
	std::set<std::string> opset;
	TestSet * tests = ((CTest &)tspace.get_project()).malloc_test_set();
	
	/* initialize oplist */
	auto beg = data.get_subsuming_operators().begin();
	auto end = data.get_subsuming_operators().end();
	while (beg != end) {
		oplist.push_back(*(beg++));
	}

	/* select bits */
	BitSeq bits(oplist.size());
	unsigned long long times = std::pow(2, oplist.size());
	for (unsigned long long k = 0; k < times; k++) {
		/* generate the composision of operators */
		this->select_operators(oplist, bits, opset);

		/* output its coverage and dominator score */
		this->gen_cov_score(data, opset, out, *tests);

		bits.increase();	/* roll to next composition */
	}

	/* return */ 
	((CTest &)tspace.get_project()).delete_test_set(tests);
}
void SuOprtWriter::write_op_scores(SOperatorSet & data, std::ostream & out) {
	std::vector<std::string> oplist;
	std::set<std::string> opset;
	TestSet * tests = ((CTest &)tspace.get_project()).malloc_test_set();

	/* initialize oplist */
	auto beg = data.get_subsuming_operators().begin();
	auto end = data.get_subsuming_operators().end();
	while (beg != end) {
		oplist.push_back(*(beg++));
	}

	out << "operator\tcoverage\tscore\n";
	for (int i = 0; i < oplist.size(); i++) {
		const std::string & op = oplist[i];
		opset.clear(); opset.insert(op);

		out << op << "\t";
		/* output its coverage and dominator score */
		this->gen_cov_score(data, opset, out, *tests);
	}

	/* return */
	((CTest &)tspace.get_project()).delete_test_set(tests);
}

double TestMachine::evaluate(const TestSet & tests) {
	/* initialization */
	const std::set<MuCluster *> & dom_mutants
		= context.get_clusters().get_subsuming();
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

	/* return */ return ((double)K) / ((double)M);
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
		if (context.get_operators().has_operator(op))
			operators.insert(op);
	}

	/* get the subsuming clusters covered by operators */
	std::set<MuCluster *> requirements;
	const MSGraph & graph = context.get_clusters().get_graph();
	for (unsigned i = 0; i < graph.size(); i++) {
		MuCluster * cluster = &(graph.get_cluster(i));
		const std::set<std::string> & ops =
			context.get_mappings().get_operators_of(*cluster);
		auto obeg = ops.begin(), oend = ops.end();
		while (obeg != oend) {
			if (operators.count(*obeg) > 0) {
				requirements.insert(cluster);
			}
			obeg++;
		}
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
/* build up mutant subsumption graph */
static void build_subsumption_graph(MSGraph & graph, 
	ScoreProducer & producer, ScoreConsumer & consumer) {
	/* declarations */
	ScoreVector * vec;
	MSGBuilder builder;

	/* compute subsumption graph */
	builder.open(graph);
	while ((vec = producer.produce()) != nullptr) {
		builder.add(vec->get_mutant(), vec->get_vector());
		consumer.consume(vec);
	}
	builder.link(MSGLinker::randomly);
	builder.close();

	/* end */ return;
}
/* build up mutant subsumption graph with specified filter */
static void build_subsumption_graph(MSGraph & graph,
	ScoreProducer & producer, ScoreConsumer & consumer,
	MutantSpace & mspace, bool(*filter)(Mutant &)) {
	/* declarations */
	ScoreVector * vec;
	MSGBuilder builder;

	/* add into the graph as node */
	builder.open(graph);
	while ((vec = producer.produce()) != nullptr) {
		/* get mutant for evaluation */
		Mutant::ID mid = vec->get_mutant();
		Mutant & mutant = mspace.get_mutant(mid);

		/* add the mutant when available */
		if (filter(mutant)) 
			builder.add(mid, vec->get_vector());

		/* delete vector */ consumer.consume(vec);
	}
	builder.link(); builder.close();

	/* return */ return;
}

/* filter methods */
bool is_trap_mutant(Mutant & mutant) {
	const std::string & op = mutant.get_operator();
	if (op == "u-STRP" || op == "u-STRI" || op == "u-VDTR")
		return true;
	else return false;
}

/* test */
int main() {
	// input-arguments
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "triangle";
	TestType ttype = TestType::general;

	// get root file and analysis dir 
	File & root = *(new File(prefix + prname));

	// create code-project, mutant-project, test-project
	CProgram & program = *(new CProgram(root));
	CTest & ctest = *(new CTest(ttype, root, program.get_exec()));
	CMutant & cmutant = *(new CMutant(root, program.get_source()));
	CScore & cscore = *(new CScore(root, cmutant, ctest));

	// load mutants and tests
	load_tests_mutants(ctest, cmutant);

	// for alpha-list select
	std::vector<double> alphas;
	alphas.push_back(0.60); 
	alphas.push_back(0.80); 
	alphas.push_back(0.90);

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

		// compute the subsumption graph
		MSGraph graph(mspace);
		build_subsumption_graph(graph, producer, consumer);
		std::cout << "Subsumption Graph: constructed finished...\n";

		// output graph-file 
		MSGraphPrinter printer; File * target = nullptr;
		const std::vector<File *> & files = root.list_files();
		for (int i = 0; i < files.size(); i++) {
			File * child = files[i];
			if (child->get_local_name() == "analysis") {
				target = child; break;
			}
		}

		/* write information */
		if (target == nullptr) {
			CError error(CErrorType::Runtime, 
				"main()", "Invalid project: no analysis directory");
			CErrorConsumer::consume(error);
		}
		else {
			printer.open(*target);
			printer.write(graph);
			printer.close();
			std::cout << "MSG-print finished...";
		}

		// relevant object for MSG 
		MuClusterSet clusters(graph);
		OpMutantMap operators(mspace);
		OpClusterMap mappings(graph);
		std::cout << "Environment building: finished...\n";

		// extract subsuming operators
		SOperatorSet sopset(clusters, mappings, operators);
		std::cout << "Subsuming operators: extracted...\n";

		// output 
		SuOprtWriter writer(ctest.get_space());
		writer.open(root);
		writer.write(sopset);
		writer.eval_operators(sopset, alphas);
		writer.close();
		std::cout << "Output: finished...\n";

		// end this file
		std::cout << "End file: \"" << cfile.get_file().get_path() << "\"\n";
	}

	// delete memory
	delete &cscore;
	delete &cmutant; delete &ctest;
	delete &program; delete &root;

	/* exit */ std::cout << "\nPress any key to exit...\n"; getchar(); exit(0);
}
