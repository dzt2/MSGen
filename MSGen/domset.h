#pragma once

/*
	-file : domset.h
	-purp : to define interface to compute dominator mutants from DMSG
	-arth : Lin Huan
	-date : Nov 2nd, 2017
	-clas :
		class ScoreMatrix
		class MutSet
		class DomSetBuilder
*/

#include "sgraph.h"

// class declarations
class ScoreMatrix;
class MutSet;
class DomSetBuilder;

class DomSetBuilder_Greedy;
class DomSetBuilder_Blocks;

/* score matrix */
class ScoreMatrix {
public:
	/* create a matrix based on the number of mutants and tests as inputs */
	ScoreMatrix(MutantSpace & ms, TestSpace & ts) : mspace(ms), tspace(ts) {
		matrix = new BitSeq(ms.number_of_mutants() * ts.number_of_tests());
		record = new BitSeq(ms.number_of_mutants() * ts.number_of_tests());
		matrix->clear_bytes(); record->clear_bytes();
	}
	/* deconstructor */
	~ScoreMatrix() { delete matrix; }
	/* add score vector into the matrix on its specified index */
	void add_score_vectors(ScoreProducer &, ScoreConsumer &);

	/* get the mutant space */
	inline MutantSpace & get_mutant_space() const { return mspace; }
	/* get the test space */
	inline TestSpace & get_test_space() const { return tspace; }
	/* get the index of result in matrix specified by mutant and test */
	inline BitSeq::size_t get_index_of(Mutant::ID mid, TestCase::ID tid) const {
		return mid * tspace.number_of_tests() + tid;
	}
	/* get the result of mutant tested against specified one */
	inline bit get_result(Mutant::ID mid, TestCase::ID tid) const {
		record->set_bit(mid * tspace.number_of_tests() + tid, BIT_1);
		return matrix->get_bit(mid * tspace.number_of_tests() + tid);
	}

	/* get the number of testings during execution */
	inline size_t number_of_testings() const {
		size_t count = 0;
		BitSeq::size_t n = record->bit_number();
		for (BitSeq::size_t k = 0; k < n; k++) {
			if (record->get_bit(k)) count++;
		}
		return count;
	}
	/* clear all the tested records */
	inline void clear_all_testings() {
		record->clear_bytes();
	}
	/* get the total number of testings (maximum) */
	inline size_t number_of_all_testings() const {
		return record->bit_number();
	}
	/* get number of equivalents */
	inline size_t get_equivalents() const { return equivalents; }

private:
	MutantSpace & mspace;
	TestSpace & tspace;
	BitSeq * matrix;
	BitSeq * record;
	size_t equivalents;
};
/* set for mutant records */
class MutSet {
public:
	/* create an empty set of mutants */
	MutSet(MutantSpace & space) : mspace(space), mutants() {}
	/* deconstructor */
	~MutSet() { clear_mutants(); }

	/* get the mutant space for this set */
	inline MutantSpace & get_space() const { return mspace; }
	/* get the number of mutants in the set */
	inline size_t size() const { return mutants.size(); }
	/* add this mutant to this set */
	inline void add_mutant(Mutant::ID mid) { mutants.insert(mid); }
	/* get the set of mutants in the set */
	inline const std::set<Mutant::ID> & get_mutants() const { return mutants; }
	/* clear all the mutants in this set */
	inline void clear_mutants() { mutants.clear(); }

private:
	MutantSpace & mspace;
	std::set<Mutant::ID> mutants;
};
/* builder for dominator set */
class DomSetBuilder {
protected:
	/* construct a builder */
	DomSetBuilder() : matrix(nullptr) {}
	/* deconstructor */
	virtual ~DomSetBuilder() { close(); }

	/* compute the dominator set from score matrix */
	virtual void compute(MutSet & ans) {
		CError error(CErrorType::InvalidArguments, 
			"DomSetBuilder::compute", 
			"Invalid access: unimplemented method...");
		CErrorConsumer::consume(error); 
		exit(CErrorType::InvalidArguments);
	}

public:
	/* open to receive the inputs */
	void open(ScoreMatrix & mat) { 
		matrix = &mat; 
	}
	/* build up the dominator set in ans based on matrix */
	void build(MutSet & ans) {
		if (matrix == nullptr) {
			CError error(CErrorType::InvalidArguments,
				"DomSetBuilder::build",
				"Invalid access: not openned");
			CErrorConsumer::consume(error);
			exit(CErrorType::InvalidArguments);
		}
		else {
			matrix->clear_all_testings();
			ans.clear_mutants(); compute(ans);
		}
	}
	/* close the builder */
	void close() { matrix = nullptr; elist.clear(); }

	/* get the number of eliminated mutants each loop */
	inline const std::vector<size_t> & get_eliminates() const { return elist; }

protected:
	/* score function as basis for computation */
	ScoreMatrix * matrix;
	/* list to record the number of mutants eliminated by each loop */
	std::vector<size_t> elist;
};

/* to determine dominator set based on classical greedy algorithm */
class DomSetBuilder_Greedy : public DomSetBuilder {
public:
	/* create a builder based on greedy algorithm */
	DomSetBuilder_Greedy() : DomSetBuilder(), 
		compare_counts(0), select_times(0) {}
	/* deconstructor */
	~DomSetBuilder_Greedy() { }

	/* get the number of comparisons */
	inline size_t get_comparisons() const { return compare_counts; }
	/* get the times to select mutants */
	inline size_t get_selections() const { return select_times; }

protected:
	/* compute the dominator set based on classical algorithm */
	void compute(MutSet & ans);

	/* get the score set for each mutant */
	void derive_score_set(Mutant::ID, std::set<TestCase::ID> &);
	/* whether the mutant is killed by all the tests in set */
	bool is_killed_by_all(Mutant::ID, const std::set<TestCase::ID> &);
	/* eliminate those subsummed by the mi from set */
	void erase_subsummeds(Mutant::ID, std::set<Mutant::ID> &);
	/* get the next unvisited mutant from M, if all are visited, return false */
	bool get_next_mutants(const std::set<Mutant::ID> & M,
		const std::set<Mutant::ID> & records, Mutant::ID & next);

private:
	/* number to compare two mutants to determine subsumption */
	size_t compare_counts;
	/* times to select mutant */
	size_t select_times;

};
/* to determin dominator set based on coverage approach */
class DomSetBuilder_Blocks : public DomSetBuilder {
public:
	/* construct the builder based on coverage-approach */
	DomSetBuilder_Blocks(MS_Graph & csg) :
		DomSetBuilder(), graph(csg),
		cover_sets(), score_sets() {}
	/* destructor */
	~DomSetBuilder_Blocks() { clear_cache(); }

	/* get the number of comparisons */
	inline size_t get_comparisons() const { return compare_counts; }
	/* get the times to select mutants */
	inline size_t get_selections() const { return select_times; }

private:
	/* graph for block dominance */
	MS_Graph & graph;

protected:
	/* compute the set of mutants */
	void compute(MutSet & ans);

	// inner-block-methods
	/* collect the set of mutants seeded within the block from M */
	void collect_mutant_of(const std::set<Mutant::ID> &, MSG_Node &, std::set<Mutant::ID> &);
	/* collect the set of tests that cover the block */
	const std::set<TestCase::ID> & derive_cover_set(MSG_Node &);
	/* get the score set of a mutant based on its coverage set */
	const std::set<TestCase::ID> & derive_score_set(const std::set<TestCase::ID> &, Mutant::ID);
	/* whether the mutant is killed by each test in the set */
	bool is_killed_by_all(Mutant::ID, const std::set<TestCase::ID> &);
	/* eliminate those subsummed by the mi from M based on its score set (mid is left behind) */
	void erase_subsummeds(Mutant::ID, std::set<Mutant::ID> &, const std::set<TestCase::ID> &);
	/* get the next unvisited mutant from M, if all are visited, return false */
	bool get_next_mutants(const std::set<Mutant::ID> & M,
		const std::set<Mutant::ID> & records, Mutant::ID & next);
	/* compute the dominator mutants within one block */
	void compute_inner_block(std::set<Mutant::ID> &, MSG_Node &, std::set<Mutant::ID> &);

	// inter-block-methods
	/* get the blocks that can be covered by all the tests of score set, in the context of given blocks */
	void get_feasible_domain(const std::set<TestCase::ID> & scoreset, std::set<MSG_Node *> & domain);
	/* get the mutants within the set of blocks from the current mutant set */
	void collect_mutants_of(const std::set<Mutant::ID> & M,
		const std::set<MSG_Node *> & blocks, std::set<Mutant::ID> & mutants);
	/*
		For each mutant in dominator set of block, do the followng job;
		1) compute the feasible domain by mi's score set;
		2) compute the mutants among the feasible domain;
		3) erase those subsummed by mi in mutants of domain;
		4) put the remainders back to the set M;
	*/
	void compute_inter_block(std::set<Mutant::ID> & M, MSG_Node & block, const std::set<Mutant::ID> domset);

	/* update the mutants by each block loop:
		1) compute the dominator set in specified block, said Di;
		2) update the mutants by inter-block relations;
		3) reput the mutants in Di back to the set M.
	*/
	void update_by_block(std::set<Mutant::ID> & M, MSG_Node & block);
	/* sort the sequence of blocks by their dominance */
	void sort_blocks_by(MS_Graph & csg, std::vector<MSG_Node *> & list);
	/* get the node with minimal number of mutants in set */
	MSG_Node * get_min_mutants(const std::set<MSG_Node *> & nodes);

	// resource-methods
	void clear_cache();
private:
	size_t compare_counts, select_times;
	std::map<Mutant::ID, std::set<TestCase::ID> *> score_sets;
	std::map<MSG_Node *, std::set<TestCase::ID> *> cover_sets;
};