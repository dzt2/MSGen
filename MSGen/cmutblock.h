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
class CMutantBlockBuilder;
class CMutantBlockMerger;

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
	/* insert mutant and modify local-graph */
	friend class CMutantBlockBuilder;

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
	/* clear all the blocks in the set */
	void clear();

private:
	/* mutant space */
	MutantSpace & mspace;
	/* set of mutant blocks */
	std::set<CMutantBlock *> blocks;
	/* index from mutant to block */
	std::map<Mutant::ID, CMutantBlock *> index;
};
/* builder for mutations blocks */
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
