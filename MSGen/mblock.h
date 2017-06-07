#pragma once

/*
	-File : mblock.h
	-Arth : Lin Huan
	-Date : Jun 7th 2017
	-Purp : to define model for mutation blocks and their builders
	-Clas :
		[1] MutBlock
		[2] MutBridge
		[3] MutBlockGraph

		[4] MutBlockBuilder
		[5] LocalMSGBuilder
		[6] MutBridgeLinker
			_BridgeComputer
*/

#include "cscore.h"
#include "msgraph.h"

class MutBlock;
class BlockPort;
class MutBlockGraph;

class MutBlockBuilder;
class LocalMSGBuilder;
class MutBlockConnect;
class MutPortComputer;

/* mutation block */
class MutBlock {
public:
	/* type for id of block */
	typedef unsigned int ID;

	/* get the id of the block */
	ID get_block_id() const { return bid; }
	/* get the vector for coverage of this block */
	const BitSeq & get_coverage() const { return coverage; }
	/* get the set of mutants in this block */
	const std::set<Mutant::ID> & get_mutants() const { return mutants; }
	/* get the local subsumption graph */
	MSGraph & get_local_graph() { return local_graph; }

	/* whether specified mutant is in this block */
	bool has_mutant(Mutant::ID & mid) const { return mutants.count(mid) > 0; }
	/* get the graph where ths block is defined */
	MutBlockGraph & get_graph() const { return graph; }

	/* whether the block is linked to another block? */
	bool is_linked_to(MutBlock & block) const { return ports.count(&block) > 0; }
	/* get the port the target block (if defined) */
	BlockPort & get_port_of(MutBlock &) const;
	/* get the map from target block to their port */
	const std::map<MutBlock *, BlockPort *> & get_ports() const { return ports; }

	/* create and delete, add mutants */
	friend class MutBlockGraph;
	/* to link|dislink blocks */
	friend class MutBlockConnect;

private:
	/* block id */
	ID bid;
	/* coverage of this vector */
	const BitSeq coverage;
	/* graph for mutant blocks */
	MutBlockGraph & graph;
	/* set of mutants in this block */
	std::set<Mutant::ID> mutants;
	/* local mutant subsumption graph */
	MSGraph & local_graph;
	/* ports to target blocks in this block */
	std::map<MutBlock *, BlockPort *> ports;

protected:
	/* create an empty block in specified graph */
	MutBlock(MutBlockGraph & g, ID id, const BitSeq & covvec);
	/* deconstructor */
	~MutBlock() { clear_mutants(); dislink_from_all(); }

	/* add a mutant into the block */
	void add_mutant(Mutant::ID mid) { mutants.insert(mid); }
	/* clear all mutants and the local graph */
	void clear_mutants() { mutants.clear(); local_graph.clear(); }

	/* link this block to the target block */
	void link_to_block(MutBlock &);
	/* disconnect this block from specified target */
	void dislink_from_block(MutBlock &);
	/* dislink the block from all the blocks */
	void dislink_from_all();

};
/* a port is the interface for vertices in block to subsume vertices in another block */
class BlockPort {
protected:
	/* create a port from source to target */
	BlockPort(MutBlock & src, MutBlock & trg) 
		: source(src), target(trg), valid_nodes(), edges() {}
	/* deconstructor */
	~BlockPort() { clear(); }

	/* initialize the edges and valid-nodes in port */
	void init();
	/* connect the edge from src to trg */
	void link(MSGVertex &, MSGVertex &);
	/* update the valid-nodes by those with subsumption in edges */
	void update();
	/* clear all nodes and edges between blocks (in port) */
	void clear();

public:
	/* get the source block */
	MutBlock & get_source_block() const { return source; }
	/* get the target block */
	MutBlock & get_target_block() const { return target; }
	/* get the valid nodes from source block , a node is valid when:
		1) there are vertices in another block subsumed by it 
	*/
	const std::set<MSGVertex *> & get_valid_nodes() const { return valid_nodes; }
	/* get the nodes directly subsumed, by source node, in target block */
	std::list<MuSubsume> & get_subsume_nodes(MSGVertex &) const;
	
	/* create and delete */
	friend class MutBlock;
	/* init-nodes, linke-nodes, update or clear */
	friend class MutPortComputer;

private:
	/* source block */
	MutBlock & source;
	/* target block */
	MutBlock & target;
	/* valid nodes in source block */
	std::set<MSGVertex *> valid_nodes;
	/* subsumption between vertices in different blocks */
	std::map<MSGVertex *, std::list<MuSubsume> *> edges;
};
/* graph for mutant blocks */
class MutBlockGraph {
public:
	/* create block graph for specified mutants in space */
	MutBlockGraph(MutantSpace &);
	/* deconstructor */
	~MutBlockGraph() { clear(); }

	/* get the space of mutations in the blocks */
	MutantSpace & get_mutant_space() const { return mspace; }
	/* get the number of blocks in the graph */
	size_t number_of_blocks() const { return blocks.size(); }
	/* whether the specified block of id is defined in the graph */
	bool has_block(MutBlock::ID bid) const { return bid < blocks.size(); }
	/* get the block of specified id in the graph */
	MutBlock & get_block(MutBlock::ID) const;
	/* get all the blocks in graph */
	const std::vector<MutBlock *> & get_blocks() const { return blocks; }

	/* whether there is block to the mutant */
	bool has_block_of_mutant(Mutant::ID mid) const { return index.count(mid) > 0; }
	/* get the block of the specified mutant */
	MutBlock & get_block_of_mutant(Mutant::ID) const;

	/* new block|add mutant|clear */
	friend class MutBlockBuilder;

private:
	/* mutant space */
	MutantSpace & mspace;
	/* set of all blocks in the graph */
	std::vector<MutBlock *> blocks;
	/* index from mutant to their block */
	std::map<Mutant::ID, MutBlock *> index;

protected:
	/* create a new block with specified coverage */
	MutBlock * new_block(const BitSeq &);
	/* add mutant to specified block */
	void add_mutant(MutBlock *, Mutant::ID);
	/* clear all blocks (and their mutants and edges) from graph */
	void clear();
	
};

/* builder to create block (with empty mutants) in graph */
class MutBlockBuilder {
public:
	/* create an empty builder for blocks */
	MutBlockBuilder();
	/* deconstructor */
	~MutBlockBuilder() { close(); }

	/* build blocks in graph with specified coverage in producer */
	void build_blocks(MutBlockGraph &, CoverageProducer &, CoverageConsumer &);

protected:
	/* open to the graph */
	void open(MutBlockGraph &);
	/* build block for the next coverage vector */
	void build(const CoverageVector &);
	/* close the builder */
	void close();

private:
	/* graph to be built up */
	MutBlockGraph * graph;
	/* trie for distinguishing coverage */
	BitTrieTree * trie;
};
/* builder for local MSG in each block */
class LocalMSGBuilder {
public:
	/* create an empty builder for local MSG */
	LocalMSGBuilder();
	/* deconstructor */
	~LocalMSGBuilder() { close(); }

	/* construct local MSG for graph */
	void construct_local_graph(MutBlockGraph &, ScoreProducer &, ScoreConsumer &);

protected:
	/* open next graph for constructing */
	void open(MutBlockGraph &);
	/* create local graph for next score vector */
	void construct(const ScoreVector &);
	/* close the builder */
	void close();

private:
	/* blocks in the graph to be constructed */
	MutBlockGraph * graph;
	/* to build local MSG */
	MSGBuilder builder;
};
/* To connect blocks */
class MutBlockConnect {
public:
	/* create an empty connect */
	MutBlockConnect();
	/* deconstructor */
	~MutBlockConnect() { close(); }
	
	/* link the blocks and subsumption in graph */
	void connect(MutBlockGraph &);
protected:
	/* open next graph for connecting */
	void open(MutBlockGraph &);
	/* create initial connection (and their valid nodes) between blocks (with intersection) */
	void init_connection();
	/* compute the subsumption and update its valid nodes in port */
	void build_up_port(BlockPort &);
	/* close the connector */
	void close();

private:
	/* blocks in the graph to be constructed */
	MutBlockGraph * graph;
	/* to compute nodes and subsumption for bridge */
	MutPortComputer * linker;
};
/* to compute the valid nodes and subsumption in port */
class MutPortComputer {
protected:
	MutPortComputer(BlockPort &);
	~MutPortComputer();

	// TODO...
};





