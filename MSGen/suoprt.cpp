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
	else return;
}
bool SOperatorSet::solve_by_greedy(
	double alpha, std::set<std::string> & target) {
	/* compute the constraint */
	unsigned bound = mappings.get_graph().
		get_space().number_of_mutants();

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
		else if (parent.sc_state.empty()) {
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
	else return;
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
	}

	const std::set<MuCluster *> & subsumings
		= clusters.get_subsuming();
	auto cbeg = subsumings.begin();
	auto cend = subsumings.end();
	while (cbeg != cend)
		root.sc_state.insert(*(cbeg++));

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

