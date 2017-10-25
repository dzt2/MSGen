#pragma once

/*
	-file : mclass.h
	-purp : to define classed node graph for MSG
	-arth : Lin Huan
	-date : oct 20th 2017
	-clas :
		typedef MType
		MS_C_Node
		MS_C_Graph
		MS_C_Build
*/

#include "sgraph.h"

class MS_C_Node;
class MS_C_Graph;
class MS_C_Build;

/* to represent the type of mutants */
typedef std::string MType;
/* node in classified MSG */
class MS_C_Node {
public:
	typedef unsigned int ID;

	/* get the graph where the node is defined */
	MS_C_Graph & get_graph() const { return graph; }
	/* get the node's id in the graph */
	ID get_node_id() const { return id; }
	/* the type to which the mutants of this node belong */
	const MType & get_type() const { return type; }

	/* mutants of this node */
	size_t number_of_mutants() const { return mutants; }
	/* get the id to the original MSG-Node */
	long get_source() const { return source; }

	/* output degree */
	size_t get_ou_degree() const { return nexts.size(); }
	/* input degree */
	size_t get_in_degree() const { return prevs.size(); }
	/* whether this node linked to target */
	bool is_linked_to(MS_C_Node & y) const { return nexts.count(&y) > 0; }
	/* nodes to this node points */
	const std::set<MS_C_Node *> & get_next_nodes() const { return nexts; }
	/* nodes from this one pointed */
	const std::set<MS_C_Node *> & get_prev_nodes() const { return prevs; }

	/* new | delete */
	friend class MS_C_Graph;

protected:
	/* create an unlinked node in classified graph */
	MS_C_Node(MS_C_Graph & g, ID nid, long src, const MType & typ)
		: graph(g), id(nid), type(typ), mutants(0) {}
	/* deconstructor */
	~MS_C_Node() { prevs.clear(); nexts.clear(); }

	/* link this node to target one (no reflexive and duplicated one) */
	bool link_to(MS_C_Node &);
	/* add mutant */
	bool add_mutant(Mutant::ID mid) { mutants++; return true; }

private:
	/* graph where this node is defined */
	MS_C_Graph & graph;
	/* id of the node */
	ID id; 
	/* id to the original MSG_Node */
	long source;
	/* type of the mutants in this node */
	MType type;
	/* nodes from this one pointed */
	std::set<MS_C_Node *> prevs;
	/* nodes to this one points */
	std::set<MS_C_Node *> nexts;
	/* number of mutants */
	size_t mutants;
};
/* graph for classified nodes */
class MS_C_Graph {
public:
	/* create an empty graph */
	MS_C_Graph() : nodes(), node_map() {}
	/* deconstructor */
	~MS_C_Graph() { clear_all(); }

	/* number of nodes in the graph */
	size_t size() const { return nodes.size(); }
	/* get the node by its id */
	MS_C_Node & get_node(size_t) const;
	/* whether the mutant refers to a node in this graph */
	bool has_node_of(Mutant::ID mid) const { return node_map.count(mid) > 0; }
	/* get the node of this mutant */
	MS_C_Node & get_node_of(Mutant::ID) const;

	/* create a new node in the graph */
	MS_C_Node & new_node(const MType &, long source);
	/* add the mutant to the node in graph */
	bool add_mutant(MS_C_Node &, Mutant::ID);
	/* connect node x to y */
	bool connect(MS_C_Node &, MS_C_Node &);
	/* delete all the nodes and edges in graph */
	void clear_all();

private:
	std::vector<MS_C_Node *> nodes;
	std::map<Mutant::ID, MS_C_Node *> node_map;
};
/* classifier for MSG nodes */
class MS_C_Build {
public:
	/* constructor */
	MS_C_Build() : typelib(nullptr), node_cnodes(), cnode_node() {}
	/* deconstructor */
	~MS_C_Build() { close(); }

	/* open the builder by giving its type lib */
	void open(const std::map<Mutant::ID, MType> & lib) {
		close(); typelib = &lib;
	}
	/* build up the classified graph from original MSG */
	void classify(const MS_Graph &, MS_C_Graph &);
	/* close the builder */
	void close() { 
		typelib = nullptr; 
		cnode_node.clear();
		node_cnode.clear();

		auto beg = node_cnodes.begin();
		auto end = node_cnodes.end();
		while (beg != end)
			delete (beg++)->second;
		node_cnodes.clear();
	}

protected:
	void collect_nodes(const MS_Graph &, MS_C_Graph &);
	void connect_nodes(MS_C_Graph &);

private:
	const std::map<Mutant::ID, MType> * typelib;
	std::map<MSG_Node *, std::set<MS_C_Node *> *> node_cnodes;
	std::map<MSG_Node *, MS_C_Node *> node_cnode;
	std::map<MS_C_Node *, MSG_Node *> cnode_node;

	const MType & get_type(Mutant::ID mid) const;
};

