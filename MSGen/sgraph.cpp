#include "sgraph.h"

MSG_Port::~MSG_Port() {
	int n = edges.size(), k;
	for (k = 0; k < n; k++) 
		delete edges[k];
	edges.clear();
}
MSG_Edge & MSG_Port::get_edge(int k) const {
	if (k < 0 || k >= edges.size()) {
		CError error(CErrorType::InvalidArguments,"MSG_Port::get_edge(k)","Invalid k: " + std::to_string(k));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *edges[k];
}
bool MSG_Port::link(MSG_Node & x, MSG_Node & y) {
	MSG_Edge * edge = new MSG_Edge(x, y);
	edges.push_back(edge); return true;
}

MSG_Node::MSG_Node(MS_Graph & g, long cid, const BitSeq & svec) : 
	graph(g), id(cid), mutants(nullptr), score_vector(svec), in_port(), ou_port() {
	mutants = g.get_space().create_set();
	score_degree = svec.degree();
}
MSG_Node::~MSG_Node() {
	graph.get_space().delete_set(mutants);
}
bool MSG_Node::link_to(MSG_Node & next) {
	ou_port.link(*this, next);
	next.in_port.link(*this, next);
	return true;
}

MSG_Node & MS_Graph::get_node(long id) const {
	if (id < 0 || id >= nodes.size()) {
		CError error(CErrorType::InvalidArguments, "MS_Graph::get_node(k)", "Invalid k: " + std::to_string(id));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *nodes[id];
}
MSG_Node & MS_Graph::get_node_of(Mutant::ID mid) const {
	if (mut_node.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, "MS_Graph::get_node_of(mid)", 
			"Undefined mid: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = mut_node.find(mid); return *(iter->second);
	}
}
MSG_Node & MS_Graph::new_node(const BitSeq & svec) {
	MSG_Node * node = new MSG_Node(*this, nodes.size(), svec);
	nodes.push_back(node); return *node;
}
bool MS_Graph::add_mutant(MSG_Node & node, Mutant::ID mid) {
	if (mut_node.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments, "MS_Graph::add_mutant(node, mid)", 
			"Invalid mid: " + std::to_string(mid));
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		node.get_mutants().add_mutant(mid);
		mut_node[mid] = &node; return true;
	}
}
bool MS_Graph::clear() {
	int k, n = nodes.size();
	for (k = 0; k < n; k++) 
		delete nodes[k];
	nodes.clear(); mut_node.clear();
	return true;
}