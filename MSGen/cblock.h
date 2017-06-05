#pragma once

/*
	-File : cblock.h
	-Arth : Lin Huan
	-Date : June 5th 2017
	-Purp : to define interface for mutation block, in which all mutants with the same coverage
	-Clas :
		[1] MutBlock
		[2] MutBlockBridge
		[3] MutBlockGraph

		[4] MutBlockBuilder;
		[5] LocalMSGBuilder;
		[6] MutBlockConnect;
		[7] LocalGraphMerge;
*/

#include "cscore.h"
#include "msgraph.h"

class MutBlock;
class MutBlockBridge;
class MutBlockGraph;

class MutBlockBuilder;
class LocalMSGBuilder;
class MutBlockConnect;
class LocalGraphMerge;

/* block = {coverage; mutations; local-graph; }*/
class MutBlock {
protected:
	/* construct block for mutants with specified coverage */
	MutBlock(const BitSeq & covvec, MutantSpace & mspace)
		: coverage(covvec), mutants(), local_graph(mspace) {}
	/* deconstructor */
	~MutBlock() { clear_mutants(); }

	/* add a new mutant into the block */
	void add_mutant(Mutant::ID mid) { mutants.insert(mid); }
	/* clear all the mutants (including local graph) */
	void clear_mutants() { local_graph.clear(); mutants.clear(); }
	/* get local graph for their modification */
	MSGraph & get_graph_for_modify() { return local_graph; }

public:
	/* get the coverage of this vector */
	const BitSeq & get_coverage() const { return coverage; }
	/* get the set of mutants in block */
	const std::set<Mutant::ID> & get_mutants() const { return mutants; }
	/* get the local graph in block */
	const MSGraph & get_local_graph() const { return local_graph; }
	/* whether specified mutant is defined in the block */
	bool has_mutant(Mutant::ID mid) const { return mutants.count(mid) > 0; }

	/* to create and delete */
	friend class MutBlockGraph;
	/* to add mutant or clear */
	friend class MutBlockBuilder;

private:
	/* get the coverage of mutant block */
	const BitSeq coverage;
	/* get the set of mutations in block */
	std::set<Mutant::ID> mutants;
	/* local subsumption graph */
	MSGraph local_graph;
};
/* direct subsumption (or indistinguishable) from nodes in source block to those in target block */
class MutBlockBridge {
protected:
	/* subsumption from nodes in source block to those in target block */
	MutBlockBridge(const MutBlock & src, const MutBlock & trg)
		: source(src), target(trg), vertices(), edges() {
		init();
	}
	/* deconstructor */
	~MutBlockBridge() { clear(); }

	/* to initialize the bridge by extracting valid nodes in source block */
	void init();
	/* link (valid) node in source block to that in target */
	void links(MSGVertex &, MSGVertex &);
	/* clear the subsumption between bridge */
	void clear();
public:
	/* get the source block */
	const MutBlock & get_source_block() const { return source; }
	/* get the target block */
	const MutBlock & get_target_block() const { return target; }
	
	/* get valid vertices */
	const std::set<MSGVertex *> & get_vertices() const { return vertices; }

	/* whether the node is valid in source vertices */
	bool has_vertex(MSGVertex & src) const { return vertices.count(&src) > 0; }
	/* whether the node has edges in bridge */
	bool has_edges(MSGVertex & src) const { return edges.count(&src) > 0; }
	/* get the number of out-degree of specified node (when in vertices) */
	size_t get_degree_of(MSGVertex &) const;
	/* get the subsuming edges of the specified node (valid when it has edge) */
	const std::list<MuSubsume> & get_edges_of(MSGVertex &) const;

	/* to create and delete */
	friend class MutBlockGraph;
private:
	/* source block */
	const MutBlock & source;
	/* target block */
	const MutBlock & target;
	/* set of valid source vertices */
	std::set<MSGVertex *> vertices;
	/* edges from source vertex to the directly subsumed nodes in target block */
	std::map<MSGVertex *, std::list<MuSubsume> *> edges;
};
/* graph for mutation blocks */
class MutBlockGraph {
public:
	/* create an empty block graph in specified mutations */
	MutBlockGraph(const MutantSpace & mspace) :
		mut_space(mspace), blocks(), index(), bridges() {}
	/* deconstructor */
	~MutBlockGraph() { clear_all(); }

	/* get the mutation space where mutants in blocks are defined */
	const MutantSpace & get_mutant_space() const { return mut_space; }
	/* get the blocks in the graph */
	const std::set<MutBlock *> & get_blocks() const { return blocks; }
	/* whether specified mutant refers to some block */
	bool has_block_of(Mutant::ID mid) const { return index.count(mid) > 0; }
	/* get the block of specified mutant */
	const MutBlock & get_block_of(Mutant::ID) const;
	/* get the set of bridges between blocks */
	const std::set<MutBlockBridge *> & get_bridges() const { return bridges; }

private:
	/* mutation space where mutants in blocks are defined */
	const MutantSpace & mut_space;
	/* set of blocks in the graph */
	std::set<MutBlock *> blocks;
	/* index from mutant ID to their block in graph */
	std::map<Mutant::ID, const MutBlock *> index;
	/* set of (valid) edge between blocks (with intersection of coverage) */
	std::set<MutBlockBridge *> bridges;

protected:
	/* create a new block for specified coverage in graph */
	MutBlock * new_block(const BitSeq &);
	/* add the mutant to the corresponding block */
	void add_index(Mutant::ID, MutBlock *);
	/* create a new bridge from source to target */
	MutBlockBridge * new_bridge(MutBlock &, MutBlock &);
	/* clear all nodes and bridges between blocks */
	void clear_all();
};


