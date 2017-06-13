#pragma once

/*
	-File : mclass.h
	-Arth : Lin Huan
	-Date : Jun 10th 2017
	-Purp : to define model for mutation classification
	-Clas : 
		[1] MuFeature
		[2] MuClass
		[3] MuClassSet
		[4] MuClassifier
*/

#include "cscore.h"

class MuClass;
class MuClassSet;
class MuClassifier;
typedef void * MuFeature;

/* mutant class */
class MuClass {
public:
	/* get the set where the class is defined */
	MuClassSet & get_classes() const { return classes; }
	/* get the feature of this class */
	const MuFeature & get_feature() const { return feature; }
	/* get the mutants in this class */
	const MutantSet & get_mutants() const { return *mutants; }

	/* whether there is mutant in this class */
	bool has_mutant(Mutant::ID mid) const { return mutants->has_mutant(mid); }
	/* get the number of mutants in this class */
	size_t size() const { return mutants->number_of_mutants(); }

	/* create | delete */
	friend class MuClassSet;
	/* clear | add-mutant */
	friend class MuClassifier;
private:
	/* set where the class is defined */
	MuClassSet & classes;
	/* feature of the class */
	const MuFeature feature;
	/* set of mutations in this class */
	MutantSet * mutants;

protected:
	/* create an empty class in set with specified feature */
	MuClass(MuClassSet & cset, MuFeature ft);
	/* deconstructor */
	~MuClass();

	/* clear all mutants in this class */
	void clear_mutants() { mutants->clear(); }
	/* add mutant in this class */
	void add_mutant(Mutant::ID mid) { mutants->add_mutant(mid); }
};
/* set of mutant classes */
class MuClassSet {
public:
	/* create an empty set for classes on mutants */
	MuClassSet(const MutantSet & mutants) : base(mutants), classes() {}
	/* deconstructor */
	~MuClassSet();

	/* get the mutants where the classes are defined */
	const MutantSet & get_mutants() const { return base; }
	/* get the classes in this set */
	const std::map<MuFeature, MuClass *> & get_classes() const { return classes; }

	/* whether there is class of specified feature */
	bool has_class(MuFeature ft) const { return classes.count(ft) > 0; }
	/* get class of the specified feature */
	MuClass & get_class(MuFeature) const;

	/* create | delete | new-class */
	friend class MuClassifier;

private:
	/* base where the classes are defined */
	const MutantSet & base;
	/* set of classes in the set */
	std::map<MuFeature, MuClass *> classes;

protected:
	/* create a new class for feature */
	MuClass * new_class(MuFeature);
};
/* classifier for mutations */
class MuClassifier {
public:
	/* classify mutants and create class set */
	void classify(MuClassSet &);

protected:
	/* create classifier */
	MuClassifier() {}
	/* deconstructor */
	virtual ~MuClassifier();
	/* compute the next mutant, including its feature */
	virtual void next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & ft) {
		throw "Invalid access to virtual classifier!";
	}
};

/* classifier for mutants by their operator */
class MuClassifierByOperator : public MuClassifier {
public:
	/* classifier for mutant operator */
	MuClassifierByOperator() : MuClassifier(), operators() {}
	/* deconstructor */
	~MuClassifierByOperator() { operators.clear(); }

private:
	std::map<std::string, std::string *> operators;

protected:
	void next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & ft);
};
/* classifier for mutants by their location */
class MuClassifierByLocation : public MuClassifier {
public:
	/* classifier for mutant operator */
	MuClassifierByLocation() : MuClassifier(), locations() {}
	/* deconstructor */
	~MuClassifierByLocation() { locations.clear(); }

private:
	std::map<std::string, CodeLocation *> locations;

protected:
	/* get the next mutant and its location */
	void next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & ft);
};
/* classifier for mutants by their coverage */
class MuClassifierByCoverage : public MuClassifier {
public:
	/* classifier for mutant operator */
	MuClassifierByCoverage() : MuClassifier(), producer(nullptr), consumer(nullptr), trie(nullptr) {}
	/* deconstructor */
	~MuClassifierByCoverage() { uninstall(); }

	/* install the classifier with coverage producer */
	void install(CoverageProducer & p, CoverageConsumer & c) {
		producer = &p; consumer = &c; trie = new BitTrieTree();
	}
	/* remove the coverage producer from consideration */
	void uninstall() { 
		producer = nullptr; consumer = nullptr; 
		if (trie != nullptr) delete trie; 
		trie = nullptr;
	}

private:
	CoverageProducer * producer;
	CoverageConsumer * consumer;
	BitTrieTree * trie;

protected:
	/* get the next mutant and its location */
	void next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & ft);
};
/* classifier for mutants by their score vector */
class MuClassifierByScore : public MuClassifier {
public:
	/* classifier for mutant operator */
	MuClassifierByScore() : MuClassifier(), producer(nullptr), consumer(nullptr), trie(nullptr) {}
	/* deconstructor */
	~MuClassifierByScore() { uninstall(); }

	/* install the classifier with coverage producer */
	void install(ScoreProducer & p, ScoreConsumer & c) {
		producer = &p; consumer = &c; trie = new BitTrieTree();
	}
	/* remove the score producer from consideration */
	void uninstall() { 
		producer = nullptr; consumer = nullptr; 
		if (trie != nullptr) delete trie;
		trie = nullptr;
	}

private:
	ScoreProducer * producer;
	ScoreConsumer * consumer;
	BitTrieTree * trie;

protected:
	/* get the next mutant and its location */
	void next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & ft);
};
