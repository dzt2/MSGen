#include "mclass.h"

bool MS_C_Node::link_to(MS_C_Node & next) {
	if (nexts.count(&next) > 0) return false;
	else { 
		nexts.insert(&next); 
		next.prevs.insert(this);
		return true; 
	}
}

MS_C_Node & MS_C_Graph::get_node(size_t k) const {
	if (k >= nodes.size()) {
		CError error(CErrorType::InvalidArguments, 
			"MS_C_Node::get_node(k)", "out of index: " + std::to_string(k));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *(nodes[k]);
}
MS_C_Node & MS_C_Graph::get_node_of(Mutant::ID mid) const {
	if (node_map.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments,
			"MS_C_Node::get_node_of(mid)", "undefined mutant: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else { auto iter = node_map.find(mid); return *(iter->second); }
}
MS_C_Node & MS_C_Graph::new_node(const MType & type, long source) {
	MS_C_Node * node = new MS_C_Node(*this, nodes.size(), source, type);
	nodes.push_back(node); return *node;
}
bool MS_C_Graph::add_mutant(MS_C_Node & node, Mutant::ID mid) {
	if (node_map.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments,
			"MS_C_Node::add_mutant(mid)", "duplicated mutant: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		node.add_mutant(mid); node_map[mid] = &node; return true;
	}
}
bool MS_C_Graph::connect(MS_C_Node & x, MS_C_Node & y) {
	return x.link_to(y);
}
void MS_C_Graph::clear_all() {
	size_t i, n = nodes.size();
	for (i = 0; i < n; i++) 
		delete nodes[i];
	nodes.clear(); node_map.clear();
}

void MS_C_Build::classify(const MS_Graph & source, MS_C_Graph & graph) {
	/* initialization */
	if (typelib == nullptr) {
		CError error(CErrorType::InvalidArguments,
			"MS_C_Build::classify()", "Unopenned builder");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else graph.clear_all();

	collect_nodes(source, graph);
	connect_nodes(graph);
}
const MType & MS_C_Build::get_type(Mutant::ID mid) const {
	if (typelib->count(mid) == 0) {
		CError error(CErrorType::InvalidArguments,
			"MS_C_Build::get_type()", "Unopenned mutant: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = typelib->find(mid); return iter->second;
	}
}
void MS_C_Build::collect_nodes(const MS_Graph & source, MS_C_Graph & graph) {
	MutantSpace & space = source.get_space();
	Mutant::ID mid, n = space.number_of_mutants();

	std::map<MType, MS_C_Node *> type_nodes;
	for (mid = 0; mid < n; mid++) {
		if (source.has_node_of(mid)) {
			/* get type-node-id key for the classified node */
			MSG_Node & node = source.get_node_of(mid);
			if (node.get_score_degree() == 0) continue;	// remove equivalent mutants

			/* get the key for classified node in graph */
			long nid = node.get_node_id();
			const MType & type = get_type(mid);
			std::string key = type + "@" + std::to_string(nid);

			/* undefined then create one */
			if (type_nodes.count(key) == 0) {
				/* create a new node in classified graph */
				MS_C_Node & cnode = graph.new_node(type, nid);
				type_nodes[key] = &cnode;
				
				/* connect the classified node with original cluster */
				if (node_cnodes.count(&node) == 0) {
					node_cnodes[&node] = new std::set<MS_C_Node *>();
					node_cnode[&node] = &cnode;
				}
				auto iter = node_cnodes.find(&node);
				(iter->second)->insert(&cnode);
				cnode_node[&cnode] = &node;
			}

			/* get the node for this mutant */
			auto iter = type_nodes.find(key);
			MS_C_Node & cnode = *(iter->second);
			graph.add_mutant(cnode, mid);
		}
	}
	type_nodes.clear();
}
void MS_C_Build::connect_nodes(MS_C_Graph & graph) {
	/* compute the represents and their corresponding origins */
	auto cbeg = node_cnodes.begin(), cend = node_cnodes.end();
	while (cbeg != cend) {
		/* get the node and its cnodes in cgraph */
		const std::set<MS_C_Node *> & cnodes = *(cbeg->second);
		MSG_Node & node = *(cbeg->first); cbeg++;
		auto beg = cnodes.begin(), end = cnodes.end();
		while (beg != end) {
			MS_C_Node & x = *(*(beg++));
			for (auto next = beg; next != end; next++) {
				MS_C_Node & y = *(*(next));
				graph.connect(x, y);
				graph.connect(y, x);
			}
		}
	}

	/* connect representations together by direct subsumption */
	auto rbeg = node_cnode.begin(), rend = node_cnode.end();
	for(;rbeg != rend; rbeg++) {
		MSG_Node & xnode = *(rbeg->first);
		MS_C_Node & xcnode = *(rbeg->second);

		const MSG_Port & port = xnode.get_ou_port();
		for (int i = 0; i < port.degree(); i++) {
			MSG_Edge & edge = port.get_edge(i);
			MSG_Node & ynode = edge.get_target();
			
			auto iter = node_cnode.find(&ynode);
			MS_C_Node & ycnode = *(iter->second);

			graph.connect(xcnode, ycnode);
		}
	}
}
