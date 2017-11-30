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

class DomSetAlgorithm_Counter;

/* score matrix */
class ScoreMatrix {
public:
	/* create a matrix based on the number of mutants and tests as inputs */
	ScoreMatrix(MutantSpace & ms, TestSpace & ts) : mspace(ms), tspace(ts) {
		matrix = new BitSeq(ms.number_of_mutants() * ts.number_of_tests());
		matrix->clear_bytes(); 
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
	/* get number of equivalents */
	inline size_t get_equivalents() const { return equivalents; }
	/* whether mutant is killed by the test */
	inline bool get_result(Mutant::ID mid, TestCase::ID tid) const {
		BitSeq::size_t bias = mid * tspace.number_of_tests();
		return matrix->get_bit(bias + tid) == BIT_1;
	}

private:
	MutantSpace & mspace;
	TestSpace & tspace;
	BitSeq * matrix;
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
typedef struct {
	Mutant::ID mid;
	size_t declines;
} MutPair;
/* to count the costs in algorithm */
class DomSetAlgorithm_Counter {

public:
	/* create a recorder with M * T space */
	DomSetAlgorithm_Counter(size_t mutants, size_t tests)
		: mnum(mutants), tnum(tests),
		mutant_test(mutants * tests), state_trans(), declines() {}
	/* deconstructor */
	~DomSetAlgorithm_Counter() { init(); }

	/* clear the data records */
	void init() {
		mutant_test.clear_bytes();
		state_trans.clear();
		declines.clear();
	}
	/* record mi testing tj */
	inline void testing(Mutant::ID mid, TestCase::ID tid) {
		mutant_test.set_bit(mid * tnum + tid, BIT_1);
	}
	/* determine whether mid is equivalent */
	inline void equivalent(Mutant::ID mid) {
		if (state_trans.count(mid) == 0)
			state_trans[mid] = 0;

		auto iter = state_trans.find(mid);
		size_t num = iter->second + 1;
		state_trans[mid] = num;
	}
	/* compute subsumption from mi to mj */
	inline void compare(Mutant::ID mid, Mutant::ID mj) {
		if (state_trans.count(mid) == 0)
			state_trans[mid] = 0;

		auto iter = state_trans.find(mid);
		size_t num = iter->second + 1;
		state_trans[mid] = num;
	}
	/* record the number of Ds for next selection */
	inline void decline(Mutant::ID mid, size_t Ds_length) {
		MutPair pair; pair.mid = mid;
		pair.declines = Ds_length;
		declines.push_back(pair);
	}

	/* get #test for mutant */
	size_t get_mutant_tests(Mutant::ID mid) const {
		size_t testing = 0;
		size_t n = tnum, bias = mid * tnum;
		for (size_t i = 0; i < tnum; i++) {
			if (mutant_test.get_bit(i + bias) == BIT_1)
				testing = testing + 1;
		}
		return testing;
	}
	/* get #trans for mutant */
	size_t get_states_trans(Mutant::ID mid) const {
		if (state_trans.count(mid) == 0)
			return 0;
		else {
			auto iter = state_trans.find(mid);
			return iter->second;
		}
	}
	/* get decline-function */
	const std::vector<MutPair> & get_declines() const { return declines; }

private:
	const size_t mnum, tnum;
	BitSeq mutant_test;
	std::map<Mutant::ID, size_t> state_trans;
	std::vector<MutPair> declines;
};

/* builder for dominator set */
class DomSetBuilder {
protected:
	/* construct a builder */
	DomSetBuilder() : 
		matrix(nullptr), counter(nullptr) {}
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
		size_t mutnum = mat.get_mutant_space().number_of_mutants();
		size_t tesnum = mat.get_test_space().number_of_tests();
		counter = new DomSetAlgorithm_Counter(mutnum, tesnum);
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
			counter->init();
			ans.clear_mutants(); 
			compute(ans);
		}
	}
	/* close the builder */
	void close() {
		if (counter != nullptr) delete counter;
		counter = nullptr;
	}

	/* get the score matrix */
	inline ScoreMatrix & get_input_matrix() const { return *matrix; }
	/* get the counter for its performance */
	inline DomSetAlgorithm_Counter & get_performance_counter() const { return *counter; }

protected:
	/* score function as basis for computation */
	ScoreMatrix * matrix;
	/* counter for performance */
	DomSetAlgorithm_Counter * counter;

};

/* to determine dominator set based on classical greedy algorithm */
class DomSetBuilder_Greedy : public DomSetBuilder {
public:
	/* create a builder based on greedy algorithm */
	DomSetBuilder_Greedy() 
		: DomSetBuilder() {}
	/* deconstructor */
	~DomSetBuilder_Greedy() { }

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
	std::map<Mutant::ID, std::set<TestCase::ID> *> score_sets;
	std::map<MSG_Node *, std::set<TestCase::ID> *> cover_sets;
};



