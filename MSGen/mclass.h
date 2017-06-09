#pragma once

/*
	-File : mclass.h
	-Arth : Lin Huan
	-Date : June 9th, 2017
	-Purp : to define model for mutation classifier
	-Clas :
		[1] MutClass
		[2] MutClassGroup
		[3] MutClassifier (abs)
			-- MutClassifierByOperator
			-- MutClassifierByLocation
			-- MutClassifierByCoverage
*/

#include "cscore.h"

class MutClass;
class MutClassGroup;
class MutClassifier;
class MutClassifierByOperator;
class MutClassifierByLocation;
class MutClassifierByCoverage;

/* type for class feature */
enum MutMetaType {
	String,
	Location,
	BitSequence,
};
/* feature for mutant class */
class MutFeature {
public:
	/* create feature for string */
	MutFeature(const std::string & name) : type(String) { content = new std::string(name); }
	/* create location-feature */
	MutFeature(const CodeLocation & loc) : type(Location) { content = new CodeLocation(loc); }
	/* create bit-string feature */
	MutFeature(const BitSeq & bits) : type(BitSequence) { content = new BitSeq(bits); }
	/* deconstructor */
	~MutFeature() {
		MutFeature **list;
		switch (type) {
		case String:
			delete ((std::string *)content); break;
		case Location:
			delete ((CodeLocation *)content); break;
		case BitSequence:
			delete ((BitSeq *)content); break;
		default:
			throw "Invalid type: " + std::to_string(type); exit(2);
		}
	}

	/* get the meta type of the class */
	MutMetaType get_meta_type() const { return type; }

	/* get the string feature */
	const std::string & get_string() const;
	/* get the location feature */
	const CodeLocation & get_location() const;
	/* get the bit-string feature */
	const BitSeq & get_bit_string() const;

private:
	MutMetaType type;
	void * content;
};
/* mutation class */
class MutClass {
protected:
	/* create an empty class in group with specified feature */
	MutClass(MutClassGroup &, MutFeature &);
	/* deconstructor */
	~MutClass();

	/* add a new mutant into the class */
	void add(Mutant::ID mid) { mutants->add_mutant(mid); }
	/* clear all the mutants in the class */
	void clear() { mutants->clear(); }

public:
	/* get the group where the class is defined */
	MutClassGroup & get_group() const { return group; }
	/* get its class-description */
	const MutFeature & get_feature() const { return *feature; }
	/* get mutants of this class */
	const MutantSet & get_mutants() const { return *mutants; }
	/* whether the mutant belongs to this class */
	bool has_mutant(Mutant::ID mid) const { return mutants->has_mutant(mid); }

	/* create | delete | add | clear */
	friend class MutClassGroup;

private:
	/* group to define this class */
	MutClassGroup & group;
	/* class description, such as operator name */
	MutFeature * feature;
	/* mutants of this class */
	MutantSet * mutants;
};
/* group to define class for mutations */
class MutClassGroup {
public:
	/* create an empty class group */
	MutClassGroup(MutantSpace &);
	/* deconstructor */
	~MutClassGroup() { clear(); }

	/* get the mutant space for this group */
	MutantSpace & get_mutant_space() const { return mspace; }
	/* get the set of classes defined in the group */
	const std::map<MutFeature *, MutClass *> & get_classes() const { return classes; }
	/* whether there are mutant class of this feature */
	bool has_class(MutFeature & ft) const { return classes.count(&ft) > 0; }
	/* get the mutant class fot this feature */
	MutClass & get_class(MutFeature &) const;

	/* whether mutant refers to some class in this group */
	bool has_class_for(Mutant::ID mid) const { return index.count(mid) > 0; }
	/* get the class of specified mutant */
	MutClass & get_class_for(Mutant::ID) const;

	/* to new-class | add-mutant | clear */
	friend class MutClassifier;
private:
	/* mutant space where the classes are defined */
	MutantSpace & mspace;
	/* set of classes */
	std::map<MutFeature *, MutClass *> classes;
	/* index from mutant to their classes */
	std::map<Mutant::ID, MutClass *> index;

protected:
	/* create a new class with specified feature */
	MutClass * new_class(MutFeature &);
	/* add mutant into the class */
	void add_mutant(MutClass *, Mutant::ID);
	/* clear all classes and their mutants */
	void clear();
};

/* classifier for mutations */
class MutClassifier {
public:
	/* classify the mutants */
	void classify() { group.clear(); classify_mutants(); }

protected:
	/* group for classification */
	MutClassGroup & group;

	/* create builder for classification */
	MutClassifier(MutClassGroup & grp) : group(grp) {}
	/* deconstructor */
	virtual ~MutClassifier() {}

	/* classify mutations in the space of the group */
	virtual void classify_mutants() { throw "Invalid access to virtual class:"; }

	/* create a new class in group */
	MutClass * new_class(MutFeature & ft) { return group.new_class(ft); }
	/* add mutant into class */
	void add_mutant(MutClass * _class, Mutant::ID mid) { group.add_mutant(_class, mid); }
};
/* classifier based on operator */
class MutClassifierByOperator : public MutClassifier {
public:
	MutClassifierByOperator(MutClassGroup & grp) : MutClassifier(grp), class_map() {}
	~MutClassifierByOperator() { class_map.clear(); }
private:
	std::map<std::string, MutClass *> class_map;
protected:
	void classify_mutants();
};
/* classifier based on location */
class MutClassifierByLocation : public MutClassifier {
public:
	MutClassifierByLocation(MutClassGroup & grp) : MutClassifier(grp), class_map() {}
	~MutClassifierByLocation() { class_map.clear(); }
private:
	std::map<std::string, MutClass *> class_map;
protected:
	void classify_mutants();
};
/* classifier by coverage */
class MutClassifierByCoverage : public MutClassifier {
public:
	MutClassifierByCoverage(MutClassGroup & grp, CoverageProducer & p, CoverageConsumer & c)
		: MutClassifier(grp), trie(), producer(p), consumer(c) {}
	~MutClassifierByCoverage() {}
private:
	BitTrieTree trie;
	CoverageProducer & producer;
	CoverageConsumer & consumer;
protected:
	void classify_mutants();
};
/* classifier by score-vector */
class MutClassifierByScore : public MutClassifier {
public:
	MutClassifierByScore(MutClassGroup & grp, ScoreProducer & p, ScoreConsumer & c)
		: MutClassifier(grp), trie(), producer(p), consumer(c) {}
	~MutClassifierByScore() {}
private:
	BitTrieTree trie;
	ScoreProducer & producer;
	ScoreConsumer & consumer;
protected:
	void classify_mutants();
};
