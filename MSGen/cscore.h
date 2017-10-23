#pragma once

/*
-File : cscore.h
-Arth : Lin Huan
-Date : May 17th, 2017
-Purp : to provide interfaces to produce/consume and access to score vectors
-Clas :
[1] ScoreVector
[2] ScoreProducer
[3] ScoreConsumer
[4] ScoreFunction	{testset; mutset; template}
[5] ScoreSpace
[5] ScoreSource		{codefile; mutspace; testspace;}
[6] CScore			{CTest; CMutant}
*/

#include "bitseq.h"
#include "cmutant.h"
#include "ctrace.h"
#include "ctest.h"

class ScoreVector;
class ScoreFunction;
class ScoreSource;
class CScore;
class ScoreProducer;
class ScoreConsumer;

/* score vector */
class ScoreVector {
protected:
	ScoreVector(const ScoreFunction & func,
		Mutant::ID id, BitSeq::size_t tnum)
		: function(func), mid(id), svec(tnum), degree(0) {}
	~ScoreVector() {}

	/* to set the bit referring to this mutant as 1 */
	bool kill(TestCase::ID);
public:
	/* get the score function of the vector */
	const ScoreFunction & get_function() const { return function; }
	/* get the mutant of the score vector */
	Mutant::ID get_mutant() const { return mid; }
	/* get the vector */
	const BitSeq & get_vector() const { return svec; }
	/* get the number of tests to kill this mutant in function */
	BitSeq::size_t get_degree() const { return degree; }

	/* create and kill-set */
	friend class FileScoreProducer;
	/* delete */
	friend class ScoreConsumer;
protected:
	const ScoreFunction & function;
	Mutant::ID mid;
	BitSeq svec;
	BitSeq::size_t degree;
};
/* template to create score vector from result.txt */
class ScoreFunction {
protected:
	ScoreFunction(const ScoreSource &, const TestSet &, const MutantSet &);
	~ScoreFunction() { bid_tid.clear(); tid_bid.clear(); }

public:
	/* get the source of score function */
	const ScoreSource & get_source() const { return source; }
	/* get the tests set of this function */
	const TestSet & get_tests() const { return tests; }
	/* get the mutants set of this function */
	const MutantSet & get_mutants() const { return mutants; }

	/* get the test id the kth bit in score vector refers to */
	TestCase::ID get_test_id_at(BitSeq::size_t) const;
	/* get the index of bit in score vector that refers to the specified test */
	BitSeq::size_t get_index_of(TestCase::ID) const;

	/* create and delete */
	friend class ScoreSource;
private:
	/* space of the function */
	const ScoreSource & source;
	/* set of tests */
	const TestSet & tests;
	/* set of mutants */
	const MutantSet & mutants;
	/* map from index of bit-string in vector to their test id in test set */
	std::vector<TestCase::ID> bid_tid;
	/* map from test id in set to their index in bit-string of score vector */
	std::map<TestCase::ID, BitSeq::size_t> tid_bid;
};
/* source of ../score/{codefile}.txt */
class ScoreSource {
protected:
	ScoreSource(const CScore &, const MutantSpace &, const File &);
	~ScoreSource();

public:
	/* get the project of this source */
	const CScore & get_project() const { return project; }
	/* get the code file referring to this source */
	const CodeFile & get_codefile() const { return codefile; }
	/* get mutant space of the source */
	const MutantSpace & get_mutant_space() const { return mspace; }
	/* get test space of the source */
	const TestSpace & get_test_space() const { return tspace; }
	/* get ../score/{codefile}.txt */
	const File & get_result_file() const { return rfile; }

	/* create a new (template) function in the score source */
	ScoreFunction * create_function(const TestSet &, const MutantSet &);
	/* whether the score function is valid for access */
	bool is_accessible(ScoreFunction *) const;
	/* delete a score function in this space */
	bool delete_function(ScoreFunction *);

	/* create and delete */
	friend class CScore;
private:
	const CScore & project;
	const CodeFile & codefile;
	const MutantSpace & mspace;
	const TestSpace & tspace;
	const File & rfile;
	/* set of score functions */
	std::set<ScoreFunction *> function_pool;
};
/* create project for score analysis */
class CScore {
public:
	/* construct score project in the Siemens Suite root directory */
	CScore(const File &, const CMutant &, const CTest &);
	/* deconstructor */
	~CScore();

	/* get ../score/ */
	const File & get_directory() const { return *dir; }
	/* get the mutant project (codefile-mutspace)+ */
	const CMutant & get_mutant_project() const { return mutproject; }
	/* get the test project */
	const CTest & get_test_project() const { return testproject; }

	/* whether there is score source corresponding to the (mutants of) code file */
	bool has_source(const CodeFile & cfile) const { return sources.count(&cfile) > 0; }
	/* get the score source corresponding to the (mutants of) code file in mut-project */
	ScoreSource & get_source(const CodeFile &) const;
	/* get the set of score sources */
	const std::map<const CodeFile *, ScoreSource *> & get_sources() const { return sources; }

private:
	const File * dir;
	const CMutant & mutproject;
	const CTest & testproject;
	std::map<const CodeFile *, ScoreSource *> sources;
};

/* to produce score vectors (abstract) */
class ScoreProducer {
protected:
	ScoreProducer() {}
	virtual ~ScoreProducer() {}

public:
	virtual ScoreVector * produce() { return nullptr; }
};
/* consumer for score vectors */
class ScoreConsumer {
public:
	/* construct a score vector consumer */
	ScoreConsumer(const ScoreFunction & func) : function(&func) {}
	~ScoreConsumer() {}

	/* consume the produced vector */
	bool consume(ScoreVector * vec) {
		if (vec == nullptr) return false;
		else {
			const ScoreFunction & func = vec->get_function();
			if (&func == function) {
				delete vec; return true;
			}
			else return false;
		}
	}
protected:
	const ScoreFunction * function;
};

/* produce score vectors from file */
class FileScoreProducer : public ScoreProducer {
public:
	/* create a producer for score function */
	FileScoreProducer(const ScoreFunction &);
	/* deconstructor */
	~FileScoreProducer() {}

	/* get score function */
	const ScoreFunction & get_function() const { return function; }
	/* produce the next vector from score function's source: ../score/xxx.txt */
	ScoreVector * produce();

protected:
	const ScoreFunction & function;
	LineReader reader;
};
/* filter | select score vector for mutants */
class ScoreFilter : public ScoreProducer {
public:
	ScoreFilter(ScoreProducer & prod, const std::set<Mutant::ID> & temp)
		: producer(prod), _template(temp) {}
	~ScoreFilter() {}

	/* produce the next vector from score function's source: ../score/xxx.txt */
	ScoreVector * produce();

private:
	ScoreProducer & producer;
	const std::set<Mutant::ID> & _template;
};

/* vector to represent the coverage for each mutant */
class CoverageVector {
protected:
	/* construct all-zeros coverage */
	CoverageVector(Mutant::ID id, BitSeq::size_t tnum)
		: mid(id), vec(tnum) {}
	/* deconstructor */
	~CoverageVector() {}

	const Mutant::ID mid;
	BitSeq vec;
public:
	/* get the id of mutant for this coverage */
	Mutant::ID get_mutant() const { return mid; }
	/* return the coverage vector */
	const BitSeq & get_coverage() const { return vec; }

	/* for create */
	friend class CoverageProducer;
	/* for delete */
	friend class CoverageConsumer;
};
/* produce score vector that represents mutant coverage */
class CoverageProducer {
public:
	CoverageProducer(const MutantSpace & ms, const FileCoverage & fs)
		: mspace(ms), filecov(fs) {
		const CodeFile & cfile1 = mspace.get_code_file();
		const CodeFile & cfile2 = filecov.get_file();
		if ((&cfile1) != (&cfile2)) {
			CError error(CErrorType::InvalidArguments, "CoverageProducer::CoverageProducer", "Inconsistent file...");
			CErrorConsumer::consume(error);
		}
		else {
			beg = 0; end = mspace.number_of_mutants();
			tnum = fs.get_tests().get_space().number_of_tests();
		}
	}
	~CoverageProducer() {}

	/* get the mutants to be produced for their coverage */
	const MutantSpace & get_mutants() const { return mspace; }
	/* get the file coverage of the source code to mutants refer */
	const FileCoverage & get_coverage() const { return filecov; }

	/* produce the next score vector corresponding to the mutant coverage */
	CoverageVector * produce();
protected:
	const MutantSpace & mspace;
	const FileCoverage & filecov;
	Mutant::ID beg, end;
	BitSeq::size_t tnum;
};
/* consumer for coverage vector */
class CoverageConsumer {
public:
	/* create consumer for coverage vector */
	CoverageConsumer() {}
	/* deconstructor */
	~CoverageConsumer() {}

	/* consume next coverage vector */
	void consume(CoverageVector * cvec) { delete cvec; }
};


