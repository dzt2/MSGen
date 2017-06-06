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
	class MutBridgeLinker;
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
	/* to construct local MSG */
	friend class LocalMSGBuilder;

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
		: source(src), target(trg), vertices(), edges() {}
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
	/* to link subsumption between clusters in different blocks */
	friend class MutBridgeLinker;
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
	MutBlockGraph(MutantSpace & mspace) :
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
	MutBlock & get_block_of(Mutant::ID) const;
	/* get the set of bridges between blocks */
	const std::set<MutBlockBridge *> & get_bridges() const { return bridges; }

	/* to new-block, add-index, and clear-all */
	friend class MutBlockBuilder;
	/* to new-clear bridges between blocks */
	friend class MutBlockConnect;
private:
	/* mutation space where mutants in blocks are defined */
	MutantSpace & mut_space;
	/* set of blocks in the graph */
	std::set<MutBlock *> blocks;
	/* index from mutant ID to their block in graph */
	std::map<Mutant::ID, MutBlock *> index;
	/* set of (valid) edge between blocks (with intersection of coverage) */
	std::set<MutBlockBridge *> bridges;

protected:
	/* create a new block for specified coverage in graph */
	MutBlock * new_block(const BitSeq &);
	/* add the mutant to the corresponding block */
	void add_index(Mutant::ID, MutBlock *);
	/* create a new bridge from source to target */
	MutBlockBridge * new_bridge(MutBlock &, MutBlock &);
	/* clear all bridges between blocks */
	void clear_bridges();
	/* clear all nodes and bridges between blocks */
	void clear_all();
};

/* to build the mutation block by coverage (with empty local MSG) */
class MutBlockBuilder {
public:
	/* construct a builder for building mutations in blocks */
	MutBlockBuilder() : bgraph(nullptr), trie(nullptr) {}
	/* deconstructor */
	~MutBlockBuilder() { close(); }

	/* create the block-graph by coverage vectors */
	void build_blocks(MutBlockGraph &, CoverageProducer &, CoverageConsumer &);

protected:
	/* open another block graph for building */
	void open(MutBlockGraph &);
	/* insert mutant into block by their coverage */
	void insert(const CoverageVector &);
	/* close the current graph for building */
	void close();

private:
	/* graph for mutation blocks to be built */
	MutBlockGraph * bgraph;
	/* binary trie tree */
	BitTrieTree * trie;
};
/* to construct local MSG for each block */
class LocalMSGBuilder {
public:
	/* create builder for local MSG in blocks */
	LocalMSGBuilder() : graph(nullptr), builders() {}
	/* deconstructor */
	~LocalMSGBuilder() { close(); }
	
	/* construct the local graph for each block */
	void construct_graphs(MutBlockGraph &, ScoreProducer &, ScoreConsumer &);

protected:
	/* open another graph for processed  */
	void open(MutBlockGraph &);
	/* create clusters and hierarchy */
	void build(const ScoreVector &);
	/* create direct subsumption between clusters */
	void end();
	/* close the gragh and clear builers */
	void close();

private:
	/* graph for processed */
	MutBlockGraph * graph;
	/* map from block to their builder that construct MSG for their mutants */
	std::map<MutBlock *, MSGBuilder *> builders;
};
/* to create direct subsumption between nodes of two blocks */
class MutBridgeLinker {
public:
	/* to create, delete, connect */
	friend class MutBlockConnect;
	/* connect the subsumption between clusters in different bridge */
	void connect(MutBlockBridge &);
protected:
	/* create a linker for connection */
	MutBridgeLinker() : bridge(nullptr) {}
	/* deconstructor */
	~MutBridgeLinker() { close(); }

	/*	1) open another bridge;
	2) sort (valid) nodes in source block;
	3) initialize solutions for leafs (leafs in another graph) */
	void open(MutBlockBridge &);
	/* link the node with directly subsuming nodes in target block */
	void link(MSGVertex *);
	/* close the bridge for blocks */
	void close();
private:

	/* the node is "leaf" in valid subgraph when:
	1) none of its children are not in the set
	*/
	bool is_leaf_vertex(MSGVertex *, const std::set<MSGVertex *> &);
	/* node is solvable, only when all its children have been solved */
	bool is_solvable(MSGVertex *, const std::set<MSGVertex *> &);
	/* compute initial DS for each leaf in source block
	1) all leafs in another MSG of block
	*/
	void get_leaf_direct_subsuming(MSGVertex *, std::set<MSGVertex *> &);

	/* initialize the DS with their children's solutions */
	void initial_direct_subsuming(MSGVertex &, std::set<MSGVertex *> &);
	/* compute DS for each valid node in source block */
	void compute_direct_subsuming(MSGVertex &, std::set<MSGVertex *> &);
	/* connect x to the nodes in set DS */
	void connect_direct_subsuming(MSGVertex &, std::set<MSGVertex *> &);

	/* whether x subsumes y */
	bool subsume(MSGVertex &, MSGVertex &);

	/* bridge to be computed */
	MutBlockBridge * bridge;
	/* sorted list for solving DS for each valid vertex in source block */
	std::queue<MSGVertex *> vex_queue;
	/* solutions for each valid node in source block */
	std::map<MSGVertex *, std::set<MSGVertex *> *> solutions;
};
/* to build up bridges with direct subsumption between blocks */
class MutBlockConnect {
public:
	/* create a builder for bridge between blocks */
	MutBlockConnect() : graph(nullptr), linker() {}
	/* deconstructor */
	~MutBlockConnect() { close(); }

	/* create bridges between graph */
	void connect(MutBlockGraph &);
protected:
	/* open another block-graph */
	void open(MutBlockGraph &);
	/* create all (valid) bridges between blocks */
	void build_all_bridges();
	/* create subsumption between all valid bridges */
	void build_all_subsume();
	/* create subsumption for specified bridge */
	void build_subsume(MutBlockBridge & bridge) { linker.connect(bridge); }
	/* close the building */
	void close();

private:
	/* graph for mutant blocks */
	MutBlockGraph * graph;
	/* linker for subsumption between blocks */
	MutBridgeLinker linker;
};
