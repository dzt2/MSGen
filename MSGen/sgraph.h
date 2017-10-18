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
#include <set>
#include <queue>

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

	const MSG_Port & get_in_port() const { return in_port; }
	const MSG_Port & get_ou_port() const { return ou_port; }

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
	/* get the number of mutants in the graph */
	size_t size_of_mutants() const { return mut_node.size(); }

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
	void connect(MSG_Node & x, MSG_Node & y) { x.link_to(y); }
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
	/* constructor */
	MSG_Build(MS_Graph & g) : graph(g), producer(nullptr), consumer(nullptr) {}

	/* graph to be constructed */
	MS_Graph & graph;
	/* producer for score vector */
	ScoreProducer * producer;
	/* consumer for score vector */
	ScoreConsumer * consumer;

	/* abstract method to build up MSG */
	virtual bool construct() {
		CError error(CErrorType::Runtime, 
			"MSG_Build::construct()", "Unimplemented method access!");
		CErrorConsumer::consume(error); exit(CErrorType::Runtime);
	}

	/* whether x subsumes y */
	bool subsume(const BitSeq & x, const BitSeq & y) { return x.subsume(y); }
	/* whether x subsumes y */
	bool subsume(MSG_Node & x, MSG_Node & y) { return x.get_score_vector().subsume(y.get_score_vector()); }

public:
	/* deconstructor */
	virtual ~MSG_Build() {}

	/* open the builder for score vector module */
	void open(ScoreProducer & producer, ScoreConsumer &);
	/* close the builder for existing vector module */
	void close();
	/* build up the graph under construction from score vectors */
	bool build();
};
/* exhaustive algorithm implement: this will not construct the graph but only evaluate its performance! */
class MSG_Build_Exhaustive : public MSG_Build {
public:
	MSG_Build_Exhaustive(MS_Graph & g) : MSG_Build(g), vectors() {}
	/* deconstructor */
	~MSG_Build_Exhaustive() { clear_vectors(); }

protected:
	/* construct the graph */
	bool construct();

	/* load the score vectors from producer */
	void load_vectors();
	/* delete score vectors by the consumer */
	void clear_vectors();

	/* get |M| */
	size_t size_of_matrix() const { return size; }
	/* get matrix of |M| * |M| for subsumption between mutants */
	bool ** get_subsumption_matrix() const { return matrix; }

private:
	bool **matrix;
	size_t size;
	std::vector<ScoreVector *> vectors;
};
/* classical algorithm implement: this will update MSG */
class MSG_Build_Classical : public MSG_Build {
public:
	/* constructor */
	MSG_Build_Classical(MS_Graph & g) : MSG_Build(g) {}
	/* deconstructor */
	~MSG_Build_Classical() { }

protected:
	/* construct the graph */
	bool construct();

	/* build up the sub-graph for given score vector */
	bool build_subgraph(const std::set<ScoreVector *> &);
	/* combine the subgraph based on their roots and leafs */
	bool combine(const std::set<MSG_Node *> & ALeafs,
		const std::set<MSG_Node *> & BRoots, std::map<MSG_Node *, std::set<MSG_Node *> *> & ans);

private:
	/* classify the nodes from original into four groups */
	bool classify_vectors(const std::set<ScoreVector *> &,
		std::set<ScoreVector *> &, std::set<ScoreVector *> &,
		std::set<ScoreVector *> &, std::set<ScoreVector *> &);
	/* build up the node in graph for the score vectors (indistinguishable) */
	MSG_Node & build_up_node(const std::set<ScoreVector *> &);
	/* derive the roots and leafs among given score vector in graph */
	bool derive_roots_leafs(const std::set<ScoreVector *> &,
		std::set<MSG_Node *> &, std::set<MSG_Node *> &);

	/* whether the next node can be accessed given its visited records */
	bool topdown_ready(MSG_Node & next, const std::set<MSG_Node *> & visits);
	/* derive nodes directly subsumed by source in graph without considering its self-subsumption or by its existing children */
	bool derive_direct_subsumed(MSG_Node & source, const std::set<MSG_Node *> &troots, std::set<MSG_Node *> & DS);
	/* update the directly subsumed ones to the true answers */
	bool update_direct_subsumed(MSG_Node & source, std::set<MSG_Node *> & DS);
	/* delete those in DS that are subsumed by another one in DS */
	bool purify_direct_subsumed(std::set<MSG_Node *> & DS);
};
/* fast algorithm implement: this will update MSG */
class MSG_Build_Fast : public MSG_Build {
public:
	/* constructor */
	MSG_Build_Fast(MS_Graph & g) : MSG_Build(g) {}
	/* deconstructor */
	~MSG_Build_Fast() { }

protected:
	/* construct the graph */
	bool construct();

	/* clustering indistinguishable mutants and create nodes */
	bool clustering();
	/* ranking clusters based on their degree */
	bool rankByDegree(std::vector<std::set<MSG_Node *> *> & H);
	/* compute the nodes directly subsumed by x in base and put them into DS */
	bool directSubsumed(MSG_Node & x, const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & DS);

private:
	bool directSubsumed_topdown(MSG_Node & x, const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & DS);
	bool directSubsumed_downtop(MSG_Node & x, const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & DS);
	bool directSubsumed_randomy(MSG_Node & x, const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & DS);
	bool purify_subsumeds(std::set<MSG_Node *> & DS);

	bool derive_roots(const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & roots);
	bool derive_leafs(const std::set<MSG_Node *> & VS, std::set<MSG_Node *> & leafs);
	
	bool available_topdown(MSG_Node &y, const std::set<MSG_Node *> & VS);
	bool available_downtop(MSG_Node &y, const std::set<MSG_Node *> & VS);

	/* remove nodes that subsume y from VS */
	bool erase_subsuming(MSG_Node & y, std::set<MSG_Node *> & VS);
	/* remove nodes that are subsumed by y from VS */
	bool erase_subsumeds(MSG_Node & y, std::set<MSG_Node *> & VS);

};