#pragma once

/*
-File : msgraph.h
-Arth : Lin Huan
-Date : May 19th 2017
-Purp : to provide interfaces to access model of mutant subsumption graph
-Clas :
[1] MuCluster
MuFeature
[2] MSGVertex
[3] MSGSubsume
[4] MSGVexIndex
[5] MSGraph
[6] MSGBuilder
*/

#include "cscore.h"
#include <queue>

class MuFeature;
class MuSubsume;
class MSGVexPort;
class MSGVertex;
class MSGraph;
class MuHierarchy;
class MSGBuilder;

/* feature of mutants contains score vector and degree */
class MuFeature {
protected:
	/* construct the feature of mutant */
	MuFeature(const BitSeq & v, BitSeq::size_t d) : vector(v), degree(d) {}
	/* deconstructor */
	~MuFeature() {}
public:
	/* get the score vector of mutant */
	const BitSeq & get_vector() const { return vector; }
	/* get the mutant degree in vector */
	BitSeq::size_t get_degree() const { return degree; }

	/* create and delete */
	friend class MSGVertex;
private:
	/* score vector */
	BitSeq vector;
	/* mutant degree */
	BitSeq::size_t degree;
};
/* (direct) subsumption between mutant */
class MuSubsume {
public:
	/* create an edge from x to y*/
	MuSubsume(const MSGVertex & s, const MSGVertex & t) : source(s), target(t) {}
	/* to copy a new edge */
	MuSubsume(const MuSubsume & edge) : source(edge.source), target(edge.target) {}
	/* deconstructor */
	~MuSubsume() {}
	/* get the node from this edge points */
	const MSGVertex & get_source() const { return source; }
	/* get the node to this edge points */
	const MSGVertex & get_target() const { return target; }

	/* create and delete */
	friend class MSGVexPort;
private:
	/* node that subsumes */
	const MSGVertex & source;
	/* node that subsumed by another */
	const MSGVertex & target;
};
/* output | input port of vertex */
class MSGVexPort {
protected:
	/* create an empty port at specified node */
	MSGVexPort(const MSGVertex & vex) : node(vex), edges() {}
	/* deconstructor */
	~MSGVexPort() { edges.clear(); }

	/* create an edge from node to the target */
	void link_to(const MSGVertex & y) {
		MuSubsume edge(node, y);
		edges.push_back(edge);
	}
	/* create an edge from source to the node */
	void link_by(const MSGVertex & x) {
		MuSubsume edge(x, node);
		edges.push_back(edge);
	}
public:
	/* get the node of this port */
	const MSGVertex & get_node() const { return node; }
	/* get the edges in this port */
	const std::list<MuSubsume> & get_edges() const { return edges; }
	/* get the degree of this port */
	size_t get_degree() const { return edges.size(); }

	/* create, delete and link */
	friend class MSGVertex;
private:
	/* node of this port */
	const MSGVertex & node;
	/* edges from this port points */
	std::list<MuSubsume> edges;
};
/* vertex of mutant subsumption graph */
class MSGVertex {
public:
	typedef unsigned int ID;

	/* get the graph where the vertex is defined */
	const MSGraph & get_graph() const { return graph; }
	/* get the id of node in graph */
	ID get_id() const { return vid; }
	/* get the cluster of this node */
	const MutantSet & get_cluster() const { return *cluster; }
	/* get the feature of this node, null when non-featured */
	const MuFeature * get_feature() const { return feature; }
	/* get the input degree */
	size_t get_in_degree() const { return in_port.get_degree(); }
	/* get the output degree */
	size_t get_ou_degree() const { return ou_port.get_degree(); }
	/* get the input port of the vertex */
	const MSGVexPort & get_in_port() const { return in_port; }
	/* get the output port of the vertex */
	const MSGVexPort & get_ou_port() const { return ou_port; }

	/* create, delete, set-id, set-feature, link ports, add mutants */
	friend class MSGraph;

private:
	/* graph where the node is defined */
	const MSGraph & graph;
	/* id of the node in graph */
	ID vid;
	/* mutant cluster */
	MutantSet * cluster;
	/* mutant feature */
	MuFeature * feature;
	/* in-port */
	MSGVexPort in_port;
	/* ou-port */
	MSGVexPort ou_port;

protected:
	/* create a non-feature, non-id node in graph */
	MSGVertex(const MSGraph & g);
	/* deconstructor */
	~MSGVertex();

	/* set the id of node when added in graph */
	void set_id(ID id) { vid = id; }
	/* add a new mutant into the cluster */
	void add_mutant(Mutant::ID mid) {
		cluster->add_mutant(mid);
	}
	/* set the feature of this node */
	void set_feature(const BitSeq & v, BitSeq::size_t d) {
		if (feature != nullptr) delete feature;
		feature = new MuFeature(v, d);
	}
	/* link the node to y */
	void link_to(const MSGVertex & y) { ou_port.link_to(y); }
	/* link x to the node */
	void link_by(const MSGVertex & x) { in_port.link_by(x); }

};
/* mutant subsumption graph */
class MSGraph {
public:
	/* create an empty graph */
	MSGraph(MutantSpace & ms)
		: space(ms), vertices(), vpool(), roots(), leafs() {}
	/* deconstructor */
	~MSGraph() { clear(); }

	/* get the mutant space on this graph is defined */
	MutantSpace & get_space() const { return space; }
	/* get the number of nodes added in graph */
	size_t number_of_vertices() const { return vertices.size(); }
	/* get the nodes that has been added in graph */
	const std::vector<MSGVertex *> & get_vertices() const { return vertices; }
	/* whether node of specified id has been added in graph */
	bool has_vertex(MSGVertex::ID vid) const { return vid < vertices.size(); }
	/* get the node (added) in graph */
	MSGVertex & get_vertex(MSGVertex::ID) const;
	/* get the nodes without precendants */
	const std::set<MSGVertex *> & get_roots() const { return roots; }
	/* get the nodes without decendants */
	const std::set<MSGVertex *> & get_leafs() const { return leafs; }

	/* new-vertex, add-cluster, connect, add-vertex, clear and update roots and leafs */
	friend class MSGBuilder;
//protected:
	/* create a new (isolated) vertex in the pool */
	MSGVertex & new_vertex(const ScoreVector &);
	/* add mutant in cluster of the node */
	bool add_cluster(MSGVertex &, const ScoreVector &);
	/* create edge between two node */
	bool connect(MSGVertex &, MSGVertex &);
	/* add the node into graph */
	bool add_vertex(MSGVertex &);
	/* remove all nodes in the graph, including their mutants, and edges */
	void clear();
private:
	/* mutant space on the graph is defined */
	MutantSpace & space;
	/* vertices added in the graph */
	std::vector<MSGVertex *> vertices;
	/* vertices created by the graph */
	std::set<MSGVertex *> vpool;
	/* nodes without precendants */
protected:
	std::set<MSGVertex *> roots;
	/* nodes without decendants */
	std::set<MSGVertex *> leafs;
};

/* hierarchy for mutations clusters */
class MuHierarchy {
public:
	/* whether there is level corresponding to this degree */
	bool has_degree(BitSeq::size_t deg) const { return deg_map.count(deg) > 0; }
	/* get the level of vertices at specified degree */
	const std::vector<MSGVertex *> & get_vertices_at(BitSeq::size_t) const;
	/* get the sorted list of mutant degrees */
	const std::vector<BitSeq::size_t> & get_degrees() const { return deg_list; }
	/* get the length of hierarchy */
	size_t length() const { return deg_list.size(); }

	/* create, delete, add|clear|sort */
	friend class MSGBuilder;
protected:
	/* create an empty hierarchy */
	MuHierarchy() : deg_list(), deg_map() {}
	/* deconstructor */
	~MuHierarchy() { clear(); }

	/* add node into the hierarchy, if non-featured, return false and do nothing */
	bool add_vertex(MSGVertex &);
	/* sort the deg_list */
	void sort_degrees();
	/* clear all levels and degrees in the hierarchy */
	void clear();
private:
	/* sorted list of mutant degrees */
	std::vector<BitSeq::size_t> deg_list;
	/* map from degree to their clusters in MSG */
	std::map<BitSeq::size_t, std::vector<MSGVertex *> *> deg_map;
};
/* iterator for vertices in MSG */
class MSGIterator {
public:
	/* initialize the iterator for vertices in graph */
	virtual void open(const MSGraph &) = 0;
	/* get the next vertex in MSG (not visited yet) */
	virtual MSGVertex * get_next() = 0;
	/* erase the node from iteration */
	virtual void erase(MSGVertex *) = 0;
	/* close the iterator (clear the space) */
	virtual void close() = 0;
};
/* to iterate 'x' only when all 'y' are visited, where x subsumes y */
class MSGIter_DownTop : public MSGIterator {
public:
	MSGIter_DownTop() : visits(), window(), qset() {}
	~MSGIter_DownTop() { close(); }

	/* initialize the iterator for vertices in graph */
	void open(const MSGraph &);
	/* get the next vertex in MSG (not visited yet) */
	MSGVertex * get_next();
	/* erase the node from iteration */
	void erase(MSGVertex *);
	/* close the iterator (clear the space) */
	void close() {
		visits.clear(); window.clear();
		while (!qset.empty()) qset.pop();
		graph = nullptr;
	}

private:
	/* graph for access */
	const MSGraph * graph;
	/* set of vertices have been erased (visited) */
	std::set<MSGVertex *> visits;
	/* set of vertices going to be iterated */
	std::set<MSGVertex *> window;
	/* queue for next iteration */
	std::queue<MSGVertex *> qset;

	/* whether the given node can be visited
	1) not erased (in visits);
	2) added in graph
	3) all children have been visited
	*/
	bool accessible(MSGVertex *);
};
/* to iterate 'y' only when all 'x' are visited, where x subsumes y */
class MSGIter_TopDown : public MSGIterator {
public:
	MSGIter_TopDown() : visits(), window(), qset() {}
	~MSGIter_TopDown() { close(); }

	/* initialize the iterator for vertices in graph */
	void open(const MSGraph &);
	/* get the next vertex in MSG (not visited yet) */
	MSGVertex * get_next();
	/* erase the node from iteration */
	void erase(MSGVertex *);
	/* close the iterator (clear the space) */
	void close() {
		visits.clear(); window.clear();
		while (!qset.empty()) qset.pop();
		graph = nullptr;
	}

private:
	/* graph for access */
	const MSGraph * graph;
	/* set of vertices have been erased (visited) */
	std::set<MSGVertex *> visits;
	/* set of vertices going to be iterated */
	std::set<MSGVertex *> window;
	/* queue for next iteration */
	std::queue<MSGVertex *> qset;

	/* whether the given node can be visited
	1) not erased (in visits);
	2) added in graph
	3) all parents have been visited
	*/
	bool accessible(MSGVertex *);
};
/* the sequence by each node to be transversed is undecidable */
class MSGIter_Random : public MSGIterator {
public:
	MSGIter_Random() : graph(nullptr), nodes() {}
	~MSGIter_Random() { close(); }

	/* initialize the iterator for vertices in graph */
	void open(const MSGraph &);
	/* get the next vertex in MSG (not visited yet) */
	MSGVertex * get_next();
	/* erase the node from iteration */
	void erase(MSGVertex *);
	/* close the iterator (clear the space) */
	void close() { nodes.clear(); graph = nullptr; }

private:
	/* graph for access */
	const MSGraph * graph;
	/* nodes to be visited */
	std::set<MSGVertex *> nodes;
};

/* builder for subsumption graph */
class MSGBuilder {
public:
	/* to build subsumption graph */
	MSGBuilder() : graph(nullptr), trie(nullptr), hierarchy() {}
	/* deconstructor */
	~MSGBuilder() { close(); }

	/* build up for MSG */
	bool build(MSGraph &, ScoreProducer &, ScoreConsumer &);

	/* get the hierarchy for analysis */
	const MuHierarchy & get_hierarchy() const { return hierarchy; }

	/* clear the original graph and hierarchy */
	bool open(MSGraph &);
	/* add score vector into the builder */
	void add_score_vector(const ScoreVector &);
	/* end to add score vectors and build the graph */
	void end_score_vectors();
	
protected:
	/* clustering the mutants and create isolated vertices (not added) */
	bool clusterMutants(ScoreProducer &, ScoreConsumer &);
	/* sort the mutant clusters by their degrees in hierarchy */
	bool rankMutantByDegree();
	/* create edges between nodes and add them into graph (update roots and leafs) */
	bool build_graph();
	/* calculate the set of nodes directly subsumed by x in VS*/
	void direct_subsumed(const MSGVertex &, std::set<MSGVertex *> &);
	/* add nodes into graph and update its roots, and leafs */
	void add_vertices(const std::vector<MSGVertex *> & level);
	/* set the graph as null, other operations are invalid for access */
	bool close() {
		graph = nullptr;
		if (trie != nullptr)
			delete trie;
		return true;
	}

private:
	/* whether x subsumes y */
	bool is_subsume(const MSGVertex &, const MSGVertex &);
	/* calculate all those subsumed by x in graph */
	bool subsumed(const MSGVertex &, std::set<MSGVertex *> &);
	/* calculate all those subsuming x in graph */
	bool subsuming(const MSGVertex &, std::set<MSGVertex *> &);

	/* compute directly subsumed nodes for x, in order of down-top */
	void direct_subsumed_downtop(const MSGVertex &, std::set<MSGVertex *> &);
	/* compute directly subsumed nodes for x, in order of top-down */
	void direct_subsumed_topdown(const MSGVertex &, std::set<MSGVertex *> &);
	/* compute directly subsumed nodes for x, in order of random */
	void direct_subsumed_random(const MSGVertex &, std::set<MSGVertex *> &);

	/* subsumption graph for building */
	MSGraph * graph;
	/* trie for classifying mutants by score vector */
	BitTrieTree * trie;
	/* hierarchy of mutations */
	MuHierarchy hierarchy;

	/* down-top iterator */
	MSGIter_DownTop dt_iter;
	/* top-down iterator */
	MSGIter_TopDown td_iter;
	/* random iterator */
	MSGIter_Random rd_iter;
};