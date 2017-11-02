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
	/* get the total number of testings (maximum) */
	inline size_t number_of_all_testings() const {
		return record->bit_number();
	}

private:
	MutantSpace & mspace;
	TestSpace & tspace;
	BitSeq * matrix;
	BitSeq * record;
};
/* set for mutant records */
class MutSet {
public:
	/* create an empty set of mutants */
	MutSet(MutantSpace & space) : mspace(space), mutants() {}
	/* deconstructor */
	~MutSet() { clear_mutants(); }

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

	/* get the score set for each mutant */
	void derive_score_set(Mutant::ID, std::set<TestCase::ID> &);
	/* whether the mutant is killed by all the tests in set */
	bool is_killed_by_all(Mutant::ID, const std::set<TestCase::ID> &);

public:
	/* open to receive the inputs */
	void open(ScoreMatrix & mat) { 
		matrix = &mat; 
		compare_counts = 0;
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
			ans.clear_mutants();
			compute(ans);
		}
	}
	/* close the builder */
	void close() { matrix = nullptr; }

	/* get the number of comparisons */
	inline size_t get_comparisons() const { return compare_counts; }

protected:
	/* score function as basis for computation */
	ScoreMatrix * matrix;
	/* number to compare two mutants to determine subsumption */
	size_t compare_counts;
};

/* to determine dominator set based on classical greedy algorithm */
class DomSetBuilder_Greedy : public DomSetBuilder {
public:
	/* create a builder based on greedy algorithm */
	DomSetBuilder_Greedy() : DomSetBuilder() {}
	/* deconstructor */
	~DomSetBuilder_Greedy() { }

protected:
	/* compute the dominator set based on classical algorithm */
	void compute(MutSet & ans);

private:
	/* eliminate those subsummed by the mi from set */
	void erase_subsummeds(Mutant::ID, std::set<Mutant::ID> &);
	/* get the next unvisited mutant from M, if all are visited, return false */
	bool get_next_mutants(const std::set<Mutant::ID> & M, const std::set<Mutant::ID> & records, Mutant::ID & next);

};
/* to determin dominator set based on coverage approach */
class DomSetBuilder_Blocks : public DomSetBuilder {
public:
	DomSetBuilder_Blocks() : DomSetBuilder() {}
	~DomSetBuilder_Blocks() {}

protected:
	void compute(MutSet & ans);

private:

};