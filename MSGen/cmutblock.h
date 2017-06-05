#pragma once

/*
	-File : cmutblock.h
	-Arth : Lin Huan
	-Date : June 5th 2017
	-Purp : to generate mutant subsumption graph according to coverage optimization.
	-Clas : 
		[1] CMutantBlock
		[2] CMutantBlockSet
		[3] CMutantBlockBuilder
		[4] CMutantBlockMerger
*/

#include "ctrace.h"
#include "msgraph.h"

// class list
class CMutantBlock;
class CMutantBlockSet;
class MutantBlockBridge;
class CMutantBlockGraph;

// class for builders
class CMutantBlockBuilder;
class LocalMSGraphBuilder;
class LocalMSGraphConnect;

/* mutblock = {coverage; mutants; local-MSG;}*/
class CMutantBlock {
protected:
	/* create an empty block for mutants in space with specified coverage */
	CMutantBlock(const BitSeq & covvec, MutantSpace & mspace) 
		: coverage(covvec), mutants(), local_graph(mspace) {}
	/* deconstructor */
	~CMutantBlock() { mutants.clear(); local_graph.clear(); }
	
	/* add a new mutant in the block */
	void add_mutant(Mutant::ID mid) { mutants.insert(mid); }
	/* get graph to builder for building local MSG */
	MSGraph & get_graph_for_build() { return local_graph; }

public:
	/* get the coverage for mutants in the block */
	const BitSeq & get_coverage() const { return coverage; }
	/* get the set of mutations in the block */
	const std::set<Mutant::ID> & get_mutants() const { return mutants; }
	/* get the local subsumption graph for mutants in the block */
	const MSGraph & get_local_graph() const { return local_graph; }
	/* whether the specified mutant is in the block */
	bool has_mutant(Mutant::ID mid) const { return mutants.count(mid) > 0; }

	/* create and delete */
	friend class CMutantBlockSet;
	/* insert mutant */
	friend class CMutantBlockBuilder;
	/* build local subsumption graph in block */
	friend class LocalMSGraphBuilder;

private:
	/* coverage vector */
	const BitSeq coverage;
	/* set of mutations (with ID) */
	std::set<Mutant::ID> mutants;
	/* local mutant subsumption graph */
	MSGraph local_graph;
};
/* set of blocks for mutations */
class CMutantBlockSet {
public:
	/* create an empty set for mutation blocks */
	CMutantBlockSet(MutantSpace & space) : mspace(space), blocks() {}
	/* deconstructor */
	~CMutantBlockSet() { clear(); }

	/* get the set of blocks */
	const std::set<CMutantBlock *> & get_blocks() const { return blocks; }
	/* get the mutant space where the blocks are defined */
	MutantSpace & get_mutant_space() const { return mspace; }

	/* whether there is block for the mutant */
	bool has_block(Mutant::ID mid) const { return index.count(mid) > 0; }
	/* get the block of specified mutation */
	CMutantBlock & get_block(Mutant::ID);

	/* set of mutants in block */
	friend class CMutantBlockBuilder;
protected:
	/* create a new block in the set */
	CMutantBlock * new_block(const BitSeq &);
	/* add mutant into the block */
	void add_mutant(Mutant::ID, CMutantBlock *);
	/* initialize the bridges between nodes */
	void build_bridges();
	/* clear all the blocks in the set */
	void clear();

private:
	/* mutant space */
	MutantSpace & mspace;
	/* set of mutant blocks */
	std::set<CMutantBlock *> blocks;
	/* index from mutant to block */
	std::map<Mutant::ID, CMutantBlock *> index;
	/* set of (valid) bridges between mutant blocks */
	std::set<MutantBlockBridge *>  bridges;
};
/* bridge between local graphs in different blocks (valid cluster) */
class MutantBlockBridge {
protected:
	/* create a bridge to represents the subsumption from block-1 to block-2 */
	MutantBlockBridge(const CMutantBlock & src, const CMutantBlock & trg)
		: source(src), target(trg), source_vertices(), target_subsumes() { init(); }
	/* deconstructor */
	~MutantBlockBridge() { clear(); }

	/* initialize the bridge between block */
	void init();
	/* link the source vertex to target vertex of another graph */
	void links(MSGVertex &, MSGVertex &);
	/* clear the subsumption form nodes in source block to those of target block */
	void clear();
public:
	/* get the source block */
	const CMutantBlock & get_source_block() const { return source; }
	/* get the target block */
	const CMutantBlock & get_target_block() const { return target; }

	/* get set of source vertices (valid) to subsume nodes in another block */
	const std::set<MSGVertex *> & get_source_vertices() const { return source_vertices; }
	/* get the map from source node to their directly subsumed node in another block */
	const std::map<MSGVertex *, std::list<MuSubsume> *> & get_targets() const { return target_subsumes; }
	
	/* whether the specified node is valid to subsume nodes in target block */
	bool has_source_vertex(MSGVertex &) const;
	/* whether specified source vertex corresponds to target nodes in another block */
	bool has_target_vertices(MSGVertex &) const;
	/* get the nodes directly subsumed by the specified node in target block */
	const std::list<MuSubsume> & get_targets_of(MSGVertex &) const;
	/* the number of nodes to be directly subsumed by specified node */
	size_t get_degree_of(MSGVertex & src) const { return get_targets_of(src).size(); }

	/* to create and delete */
	friend class CMutantBlockGraph;
	/* to link nodes between blocks */
	friend class LocalMSGraphConnect;
private:
	/* source block */
	const CMutantBlock & source;
	/* target block */
	const CMutantBlock & target;
	/* vertices (valid) to subsume nodes in another block */
	std::set<MSGVertex *> source_vertices;
	/* map from source vertices to their directly subsuming nodes in target block */
	std::map<MSGVertex *, std::list<MuSubsume> *> target_subsumes;
};
/* graph between blocks (with valid connection) */
class CMutantBlockGraph {
public:
	CMutantBlockGraph();
};

/* to build mutant block (with empty local graph) */
class CMutantBlockBuilder {
public: 
	/* create a builder for block set */
	CMutantBlockBuilder() : bset(nullptr), trie(nullptr) {}
	/* deconstructor */
	~CMutantBlockBuilder() { close(); }

	/* build up the block set by coverage */
	void build_for(CMutantBlockSet &, CoverageProducer &, CoverageConsumer &);

protected:
	/* open a new block set */
	void open(CMutantBlockSet &);
	/* insert new block by coverage and mutant id */
	void insert(const CoverageVector &);
	/* close the builder for next block set */
	void close();

private:
	/* set for mutant blocks */
	CMutantBlockSet * bset;
	/* trie for searching block by coverage */
	BitTrieTree * trie;
};
/* to build up the local MSG for each block in block set */
class LocalMSGraphBuilder {
public:
	/* builder for local MSG in blocks */
	LocalMSGraphBuilder() : bset(nullptr), builders() {}
	/* deconstructor */
	~LocalMSGraphBuilder() { close(); }

	/* build up the local MSG for each block in set */
	void build_local_graph(CMutantBlockSet &, ScoreProducer &, ScoreConsumer &);

protected:
	/* open the */
	void open(CMutantBlockSet &);
	/* insert a next score vector for building local MSG */
	void insert(const ScoreVector &);
	/* close the openned block set */
	void close();

private:
	/* set of mutant blocks */
	CMutantBlockSet * bset;
	/* map from block to the builder for its local MSG */
	std::map<CMutantBlock *, MSGBuilder *> builders;
};
