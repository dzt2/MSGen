#pragma once

/*
	file : sgraph.h
	purp : to define interface for MSG and algorithms to construct it
	arth : Lin Huan
	date : oct 20th 2017
	clas :
		[1] MSG_Node
		[2] MSG_Edge
		[3] MSG_Port
		[4] MS_Graph
*/

#include "cscore.h"

class MSG_Node;
class MSG_Edge;
class MSG_Port;
class MS_Graph;
class MSG_Build;

/* edge in MSG */
class MSG_Edge {
protected:
	/* create an edge from x to y */
	MSG_Edge(MSG_Node & x, MSG_Node & y) : source(x), target(y) {}
	/* deconstructor */
	~MSG_Edge() {}

public:
	/* get the node as source in edge */
	MSG_Node & get_source() const { return source; }
	/* get the node as target in edge */
	MSG_Node & get_target() const { return target; }

	friend class MSG_Port;
private:
	MSG_Node & source;
	MSG_Node & target;
};
/* set of edges from | to a node in MSG */
class MSG_Port {
protected:
	/* create an empty collection for edges from | to a node in MSG */
	MSG_Port() : edges() {}
	/* deconstructor */
	~MSG_Port();

	/* create an edge from x to y */
	bool link(MSG_Node &, MSG_Node &);

public:
	/* number of edges in the port */
	size_t degree() const { return edges.size(); }
	/* get the kth edge in the port */
	MSG_Edge & get_edge(int k) const;

	/* create and delete */
	friend class MSG_Node;
private:
	std::vector<MSG_Edge *> edges;
};

/* node in MSG */
class MSG_Node {
protected:
	/* create an empty node with none mutants in the graph */
	MSG_Node(MS_Graph &, long, const BitSeq &);
	/* deconstructor */
	~MSG_Node();

	/* create an edge from this node to the target one */
	bool link_to(MSG_Node &);

public:
	/* get the graph where the node is defined */
	MS_Graph & get_graph() const { return graph; }
	/* get the integer id of this node in graph */
	long get_node_id() const { return id; }
	/* get the mutants belong to this node */
	MutantSet & get_mutants() const { return *mutants; }

	/* score vector for mutants in this node */
	const BitSeq & get_score_vector() const { return score_vector; }
	/* score degree for mutants in this node */
	BitSeq::size_t get_score_degree() const { return score_degree; }

	const MSG_Port get_in_port() const { return in_port; }
	const MSG_Port get_ou_port() const { return ou_port; }

	/* create | delete | link_to */
	friend class MS_Graph;

private:
	/* graph where the node is defined */
	MS_Graph & graph;
	/* node id */
	long id;
	/* mutants belong to this node */
	MutantSet * mutants;

	/* edges out from this node */
	MSG_Port ou_port;
	/* edges into this node */
	MSG_Port in_port;

	/* score vector for this node */
	const BitSeq score_vector;
	/* score degree for this node */
	BitSeq::size_t score_degree;
};
/* Mutant subsumption graph */
class MS_Graph {
public:
	/* create an empty MSG */
	MS_Graph(MutantSpace & space) : mspace(space), nodes(), mut_node() {}
	/* deconstructor */
	~MS_Graph() { clear(); }

	/* get the mutant space where MSG is defined on */
	MutantSpace & get_space() const { return mspace; }
	/* get the number of nodes in the graph */
	size_t size() const { return nodes.size(); }
	/* get the node by its id in the graph */
	MSG_Node & get_node(long) const;
	/* whether there is node for the mutant in the graph */
	bool has_node_of(Mutant::ID mid) const { return mut_node.count(mid) > 0; }
	/* get the node to which the mutant belongs */
	MSG_Node & get_node_of(Mutant::ID) const;

	/* create a new node with specified bit-string */
	MSG_Node & new_node(const BitSeq &);
	/* add mutant to the corresponding node in graph */
	bool add_mutant(MSG_Node &, Mutant::ID);
	/* connect the edge from source to target */
	bool connect(MSG_Node & x, MSG_Node & y) { x.link_to(y); }
	/* delete all the nodes and edges in the graph */
	bool clear();

private:
	MutantSpace & mspace;
	std::vector<MSG_Node *> nodes;
	std::map<Mutant::ID, MSG_Node *> mut_node;
};

/* implement algorithm to construct MSG */
class MSG_Build {
protected:
	MSG_Build(MS_Graph & g) : graph(g), producer(nullptr), consumer(nullptr) {}
	virtual ~MSG_Build() { close(); }

	MS_Graph & graph;
	ScoreProducer * producer;
	ScoreConsumer * consumer;

	virtual bool construct() {
		CError error(CErrorType::Runtime, "MSG_Build::construct()", "Unimplemented method access!");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}

public:
	void open(ScoreProducer & producer, ScoreConsumer &);
	void close();
	bool build();
};

class MSG_Build_Exhaustive : public MSG_Build {

};