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

#include "mclass.h"

// class declarations 
class MuCluster;
class MuSubsume;
class MuSubsumePort;
class MuHierarchy;
class MSGraph;
class MSGBuilder;

/* mutant cluster */
class MuCluster {
public:
	/* type for mutant cluster */
	typedef unsigned int ID;

	/* get the graph of this cluster */
	MSGraph & get_graph() const { return graph; }
	/* get the id of this graph */
	MuCluster::ID get_id() const { return cluster_id; }
	/* get the class where mutants of this cluster are maintained */
	MuClass & get_class() const { return _class; }
	/* get the score vector of this cluster */
	const BitSeq & get_score_vector() const { return score_vector; }
	/* get the score degree of this cluster */
	size_t get_score_degree() const { return score_degree; }

	/* whether there is mutant in the cluster */
	bool has_mutant(Mutant::ID mid) const { return _class.has_mutant(mid); }
	/* get the number of mutants in this cluster */
	size_t size() const { return _class.size(); }

	/* get the port of edges to this cluster */
	const MuSubsumePort & get_in_port() const { return *in_port; }
	/* get the port of edges from the cluster */
	const MuSubsumePort & get_ou_port() const { return *ou_port; }

	/* create | delete | link_to */
	friend class MSGraph;

private:
	/* graph where this cluster is defined */
	MSGraph & graph;
	/* id of this cluster */
	ID cluster_id;

	/* class where mutant is stored */
	MuClass & _class;

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
	MuCluster(MSGraph &, ID, const BitSeq &, MuClass &);
	/* deconstructor */
	~MuCluster();

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

	/* create | delete | clear | link */
	friend class MuCluster;

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
	MSGraph() : clusters(), hierarchy(), roots(), leafs(), index() {}
	/* deconstructor */
	~MSGraph() { clear(); }

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

	/* add-class | sort_hierarchy | connect | update-leafs */
	friend class MSGBuilder;

private:
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
	/* add classes in set to the MSG */
	void add(MuClassSet &);
	/* connect the edge from x to y */
	void connect(MuCluster & x, MuCluster & y) { x.link_to(y); }
	/* update the set of roots and leafs in graph */
	void update_roots_and_leafs();
};





