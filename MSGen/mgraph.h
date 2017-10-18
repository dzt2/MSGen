#pragma once

/*
	-File : mgraph.h
	-Arth : Lin Huan
	-Date : Jun 10th 2017
	-Purp : to define model for class graph (subsumption)
	-Clas :
		[1] MuCluster
		[2] MuSubsume
			MuSubsumePort
		[3] MSGraph
			MuHierarchy
			MSGIterator
		[4] MSGBuilder
*/

#include "cscore.h"
#include <queue>

// class declarations 
class MuCluster;
class MuSubsume;
class MuSubsumePort;
class MuHierarchy;
class MSGraph;
class MSGraphPrinter;

// class for building MSG
class _MSG_VSpace;
class MSGLinker;
class MSGBuilder;

/* to count the number of comparions between mutants */
extern unsigned int times;

/* mutant cluster */
class MuCluster {
public:
	/* type for mutant cluster */
	typedef unsigned int ID;

	/* get the graph of this cluster */
	MSGraph & get_graph() const { return graph; }
	/* get the id of this graph */
	MuCluster::ID get_id() const { return cluster_id; }

	/* get the mutant set of this cluster */
	const MutantSet & get_mutants() const { return *mutants; }
	/* whether there is mutant in the cluster */
	bool has_mutant(Mutant::ID mid) const { return mutants->has_mutant(mid); }
	/* get the number of mutants in this cluster */
	size_t size() const { return mutants->number_of_mutants(); }

	/* get the score vector of this cluster */
	const BitSeq & get_score_vector() const { return score_vector; }
	/* get the score degree of this cluster */
	size_t get_score_degree() const { return score_degree; }

	/* get the port of edges to this cluster */
	MuSubsumePort & get_in_port() const { return *in_port; }
	/* get the port of edges from the cluster */
	MuSubsumePort & get_ou_port() const { return *ou_port; }

	/* create | delete | link_to | add-mutant */
	friend class MSGraph;

private:
	/* graph where this cluster is defined */
	MSGraph & graph;
	/* id of this cluster */
	ID cluster_id;

	/* class where mutant is stored */
	MutantSet * mutants;

	/* score vector of the cluster */
	const BitSeq score_vector;
	/* score degree of this cluster */
	size_t score_degree;

	/* port for subsumption to this cluster */
	MuSubsumePort * in_port;
	/* port for subsumption from the cluster */
	MuSubsumePort * ou_port;

protected:
	/* create a cluster with specified cluster */
	MuCluster(MSGraph &, ID, const BitSeq &);
	/* deconstructor */
	~MuCluster();

	/* add mutant to the cluster */
	void add_mutant(Mutant::ID mid) { mutants->add_mutant(mid); }
	/* link this node to another node */
	void link_to(MuCluster &);
};
/* subsumption between cluster (direct + strict) */
class MuSubsume {
protected:
	/* create edge from source to target */
	MuSubsume(MuCluster & src, MuCluster & trg) : source(src), target(trg) {}

public:
	/* copy method */
	MuSubsume(const MuSubsume & edge) : source(edge.get_source()), target(edge.get_target()) {}
	/* deconstructor */
	~MuSubsume() {}

	/* get the source node */
	MuCluster & get_source() const { return source; }
	/* get the target node */
	MuCluster & get_target() const { return target; }

	/* create | delete */
	friend class MuSubsumePort;

private:
	/* node from which this edge points */
	MuCluster & source;
	/* node to which this edge points */
	MuCluster & target;
};
/* port to manage subsumption in cluster */
class MuSubsumePort {
protected:
	/* create an empty port */
	MuSubsumePort() : edges() {}
	/* deconstructor */
	~MuSubsumePort() { clear(); }

	/* clear all edges in the port */
	void clear() { edges.clear(); }
	/* create an edge from src to trg in the port */
	void link(MuCluster & src, MuCluster & trg) {
		MuSubsume edge(src, trg);
		edges.push_back(edge);
	}

public:
	/* get the number of edges in the port */
	size_t get_degree() const { return edges.size(); }
	/* get the list of edges in the port */
	const std::vector<MuSubsume> & get_edges() const { return edges; }

	/* create | delete | link */
	friend class MuCluster;
	/* clear */
	friend class MSGraph;

private:
	std::vector<MuSubsume> edges;
};
/* mutant hierarchy */
class MuHierarchy {
protected:
	/* create an empty hierarchy */
	MuHierarchy() : degree_list(), degree_map() {}
	/* deconstructor */
	~MuHierarchy() { clear(); }

	/* add a new cluster in the hierarchy */
	void add(MuCluster &);
	/* clear the clusters and degrees from the hierarchy */
	void clear();
	/* sort the degree list */
	void sort();

public:
	/* get the number of degress in hierarchy */
	size_t size_of_degress() const { return degree_list.size(); }
	/* get the list of mutant degrees */
	const std::vector<size_t> & get_degrees() const { return degree_list; }
	/* get the set of clusters of specified degree */
	const std::set<MuCluster *> & get_clusters_of(size_t) const;
	/* get the set of clusters of specified index in list */
	const std::set<MuCluster *> & get_clusters_at(size_t) const;

	/* create | delete | add-cluster | sort-degree */
	friend class MSGraph;

private:
	/* list of mutant degrees */
	std::vector<size_t> degree_list;
	/* map from degree to clusters */
	std::map<size_t, std::set<MuCluster *> *> degree_map;
};
/* subsumption graph */
class MSGraph {
public:
	/* create an empty graph */
	MSGraph(MutantSpace & space) : mspace(space), clusters(), hierarchy(), roots(), leafs(), index() { mutants = mspace.create_set(); }
	/* deconstructor */
	~MSGraph() { clear(); mspace.delete_set(mutants); }

	/* get the space where the graph is defined */
	MutantSpace & get_space() const { return mspace; }
	/* get the set of mutants in this graph */
	const MutantSet & get_mutants() const { return *mutants; }

	/* get the number of clusters in the grpah */
	size_t size() const { return clusters.size(); }
	/* whether there are cluster for this id in the graph */
	bool has_cluster(MuCluster::ID cid) const { return cid < clusters.size(); }
	/* get the cluster by its id in the graph */
	MuCluster & get_cluster(MuCluster::ID) const;
	/* get the mutant hierarchy */
	const MuHierarchy & get_hierarchy() const { return hierarchy; }

	/* get the roots of this graph (not subsumed) */
	const std::set<MuCluster *> & get_roots() const { return roots; }
	/* get the leafs of this graph (not subsuming) */
	const std::set<MuCluster *> & get_leafs() const { return leafs; }

	/* get the index from mutant to its cluster */
	const std::map<Mutant::ID, MuCluster *> & get_index() const { return index; }
	/* whether the mutant refers to some cluster in graph */
	bool has_cluster_of(Mutant::ID mid) const { return index.count(mid) > 0; }
	/* get the cluster of the graph */
	MuCluster & get_cluster_of(Mutant::ID) const;

	/* clear | build */
	friend class MSGBuilder;
	/* connect | update-leafs | clear-edges */
	friend class MSGLinker;

private:
	/* space where the graph is defined */
	MutantSpace & mspace;
	/* set of all mutants in the graph */
	MutantSet * mutants;

	/* set of nodes in the graph */
	std::vector<MuCluster *> clusters;
	/* hierarchy for mutant clusters */
	MuHierarchy hierarchy;

	/* set of roots in the graph */
	std::set<MuCluster *> roots;
	/* set of leafs in the graph */
	std::set<MuCluster *> leafs;

	/* map from mutant to their cluster */
	std::map<Mutant::ID, MuCluster *> index;

protected:
	/* clear all the nodes, edges and hierarchy */
	void clear();
	/* create a new cluster for bit-string in this graph */
	MuCluster * new_cluster(const BitSeq &);
	/* add mutant to the specified cluster in graph */
	void add_mutant(MuCluster &, Mutant::ID);
	/* sort the hierarchy in the graph */
	void sort() { hierarchy.sort(); }

	/* clear all the edges in the graph */
	void clear_edges();
	/* connect the edge from x to y */
	void connect(MuCluster & x, MuCluster & y) { x.link_to(y); }
	/* update the set of roots and leafs in graph */
	void update_roots_and_leafs();
};
/* print the information in graph */
class MSGraphPrinter {
public:
	/* create printer for MSG */
	MSGraphPrinter() : dir(nullptr) {}
	/* deconstructor */
	~MSGraphPrinter() { close(); }

	/* open directory where files are written */
	void open(File & file) { close(); dir = &file; }
	/* write graph information to graph.txt and mutantLib.txt in directory that has been openned */
	void write(MSGraph &);
	/* close the directory for printing */
	void close() { dir = nullptr; }

private:
	File * dir;		/* directory where graph.txt and mutantLib.txt are written */
	void replace_text(std::string &);

protected:
	/* write dir/graph.txt */
	void write_mutant_lib(MSGraph &, std::ostream &);
	/* write dir/mutantLib.txt */
	void write_mutant_graph(MSGraph &, std::ostream &);
};

/* visit space for sub-graph in MSG */
class _MSG_VSpace {
public:
	/* initialize the visit space */
	virtual void initial() {}
	/* get the next unvisited node in sub-graph */
	virtual MuCluster * next() { return nullptr; }
	/* tag all those subsuming x in sub-graph as visited */
	virtual void visit_subsuming(MuCluster *x) {}
	/* tag all those subsumed by x in sub-graph as visited */
	virtual void visit_subsumed(MuCluster *x) {}
};
/* to iterate unvisited nodes in subgraph of MSG from leafs to roots */
class _MSG_VSpace_down_top : public _MSG_VSpace {
public: 
	/* initialize the visit space */
	void initial();
	/* get the next unvisited node in sub-graph */
	MuCluster * next();
	/* tag all those subsuming x in sub-graph as visited */
	void visit_subsuming(MuCluster *x);
	/* tag all those subsumed by x in sub-graph as visited */
	void visit_subsumed(MuCluster *x) { /* do nothng! */ }

	/* create and delete */
	friend class MSGLinker;
protected:
	/* create a visit space for sub-graph from leafs to the roots */
	_MSG_VSpace_down_top(const std::set<MuCluster *> & as, 
		const std::set<MuCluster *> & lf) : adset(as), leafs(lf) {}
	/* deconstructor */
	~_MSG_VSpace_down_top();

private:
	/* set of nodes in sub-graph */
	const std::set<MuCluster *> & adset;
	/* leafs in sub-graph */
	const std::set<MuCluster *> & leafs;

	/* queue for visited node in current iteration */
	std::queue<MuCluster *> vqueue;
	/* cache for updatng vqueue in next iteration */
	std::queue<MuCluster *> vcache;
	/* set of nodes that have been (tagged as) visited */
	std::set<MuCluster *> visitset;

	/* a node is accessible when:
		1) it is in the subgraph
		2) it is not visited yet
		3) all of its children have been visited
	*/
	bool accessible(MuCluster *);
};
/* to iterate unvisited nodes in subgraph of MSG from roots to leafs */
class _MSG_VSpace_top_down : public _MSG_VSpace {
public:
	/* initialize the visit space */
	void initial();
	/* get the next unvisited node in sub-graph */
	MuCluster * next();
	/* tag all those subsuming x in sub-graph as visited */
	void visit_subsuming(MuCluster *x) { /* do nothing */ }
	/* tag all those subsumed by x in sub-graph as visited */
	void visit_subsumed(MuCluster *x);

	/* create and delete */
	friend class MSGLinker;
protected:
	/* create a visit space for sub-graph from leafs to the roots */
	_MSG_VSpace_top_down(const std::set<MuCluster *> & as,
		const std::set<MuCluster *> & rt) : adset(as), roots(rt) {}
	/* deconstructor */
	~_MSG_VSpace_top_down();

private:
	/* set of nodes in sub-graph */
	const std::set<MuCluster *> & adset;
	/* leafs in sub-graph */
	const std::set<MuCluster *> & roots;

	/* queue for visited node in current iteration */
	std::queue<MuCluster *> vqueue;
	/* cache for updatng vqueue in next iteration */
	std::queue<MuCluster *> vcache;
	/* set of nodes that have been (tagged as) visited */
	std::set<MuCluster *> visitset;

	/* a node is accessible when:
		1) it is in the subgraph
		2) it is not visited yet
		3) all of its parents have been visited
	*/
	bool accessible(MuCluster *);
};
/* to iterate unvisited nodes in subgraph of MSG in a random orders */
class _MSG_VSpace_randomly : public _MSG_VSpace {
public:
	/* initialize the visit space */
	void initial();
	/* get the next unvisited node in sub-graph */
	MuCluster * next();
	/* tag all those subsuming x in sub-graph as visited */
	void visit_subsuming(MuCluster *x);
	/* tag all those subsumed by x in sub-graph as visited */
	void visit_subsumed(MuCluster *x);

	/* create and delete */
	friend class MSGLinker;
protected:
	/* create a visit space for sub-graph from leafs to the roots */
	_MSG_VSpace_randomly(const std::set<MuCluster *> & as) : adset(as) {}
	/* deconstructor */
	~_MSG_VSpace_randomly();

private:
	/* set of nodes in sub-graph */
	const std::set<MuCluster *> & adset;
	/* set of nodes have not been visited */
	std::set<MuCluster *> visitset;
};
/* linker for direct subsumption in MSG */
class MSGLinker {
public:
	typedef char OrderOption;
	static const OrderOption down_top = 0;	/* order from leafs to roots */
	static const OrderOption top_down = 1;	/* order from roots to leafs */
	static const OrderOption randomly = 2;	/* randomly transversal */

	/* to create and connect */
	friend class MSGBuilder;
protected:
	/* create a linker to connect MSG */
	MSGLinker() : graph(nullptr), vspace(nullptr) {}
	/* deconstructor */
	~MSGLinker() { close(); }

	/* compute the direct subsumption between clusters in MSG */
	void connect(MSGraph & g) { connect(g, down_top); }
	/*	
		Compute the direct subsumption between clusters in MSG with specified orders.
		Note: this order only impacts on the efficiency, not the final results.
	*/
	void connect(MSGraph &, OrderOption);

	/* open the linker to another graph */
	void open(MSGraph &, OrderOption);
	/* compute the nodes directly subsumed by x in sub-graph */
	void compute_direct_subsumption(MuCluster &, std::set<MuCluster *> &);
	/* connect x to the nodes in DS */
	void connect_nodes(MuCluster &, const std::set<MuCluster *> &);
	/* add all the nodes in sub-graph and update its leafs | roots */
	void add_nodes_in(const std::set<MuCluster *> &);
	/* clear the sub-graph and visit-space */
	void close();

private:
	/* graph to be processed */
	MSGraph * graph;
	/* visit space for computing DS */
	_MSG_VSpace * vspace;

	/* set of nodes in sub-graph */
	std::set<MuCluster *> adset;
	/* set of leafs in sub-graph */
	std::set<MuCluster *> leafs;
	/* set of roots in sub-graph */
	std::set<MuCluster *> roots;

	/* whether x subsumes y */
	bool subsume(MuCluster &, MuCluster &);
	/* eliminate redundant mutants from DS */
	void eliminate_DS(std::set<MuCluster *> &);
};
/* to build unlinked clusters in MSG */
class MSGBuilder {
public:
	/* create a builder for MSG */
	MSGBuilder() : graph(nullptr), trie(nullptr) {}
	/* deconstructor */
	~MSGBuilder() { close(); }

	/* open another graph for building */
	void open(MSGraph &);
	/* add mutant and its score vector to the graph */
	void add(Mutant::ID, const BitSeq &);
	/* create the edges between nodes without specifying its order */
	void link();
	/* create the edges between nodes in MSG */
	void link(MSGLinker::OrderOption);
	/* close the builder and clear trie and linker */
	void close();

private:
	/* graph to be built */
	MSGraph * graph;
	/* trie tree for clustering */
	BitTrieTree * trie;
	/* to connect nodes in graph */
	MSGLinker linker;
};