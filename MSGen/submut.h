#pragma once

/*
	-File: submut.h
	-Date: July 6th 2017
	-Purp: to define interface for subsuming mutants and analysis
	-Clas:
		[1] class MutGroup
*/

#include "mgraph.h"

/* class declaration */
class MutGroup;
class MutLevel;
class MutAnalyzer;
class MutOutputer;
class TestMachine;

/* category for mutant */
typedef char MutCategory;
const MutCategory NotBelong_Category = 0;
const MutCategory Equvalent_Category = 1;
const MutCategory Subsuming_Category = 2;
const MutCategory Subsumed_Category = 3;

/* mutation group */
class MutGroup {
public:
	/* get mutant space where the mutants in group are defined */
	MutantSpace & get_space() const { return mspace; }
	/* get the mutant subsumption graph */
	const MSGraph & get_graph() const { return graph; }

	/* get the set of all mutants */
	const MutantSet & get_mutants() const { return graph.get_mutants(); }
	/* get the cluster for equivalent mutants */
	MuCluster * get_equivalents() const { return equivalents; }
	/* get the clusters for subsuming mutants */
	const std::set<MuCluster *> & get_subsumings() const { return subsumings; }

	/* get the category of the mutant in this group */
	MutCategory category_of(Mutant::ID) const;
	/* get the cluster of mutant */
	MuCluster & get_cluster_of(Mutant::ID mid) const { 
		return graph.get_cluster_of(mid); 
	}

	/* create | delete | add | classify */
	friend class MutLevel;

protected:
	/* create an empty mutant-group */
	MutGroup(MutantSpace & space) : mspace(space), graph(space),
		builder(), equivalents(nullptr), subsumings() {
		builder.open(graph);
	}
	/* deconstructor */
	~MutGroup() { clear(); builder.close(); }

	/* clear the graph, re-open builder, and clear sets */
	void clear();
	/* add mutant into the graph via builder */
	void add_mutant(Mutant::ID, const BitSeq &);
	/* classify the mutants by category */
	void classify();

private:
	/* mutant space */
	MutantSpace & mspace;
	/* subsumption graph */
	MSGraph graph;
	/* builder for MSG */
	MSGBuilder builder;

	/* cluster for equivalent mutants */
	MuCluster * equivalents;
	/* clusters for subsuming mutants */
	std::set<MuCluster *> subsumings;
};
/* used in experiment, deliver subsuming by operator -- global steps */
class MutLevel {
public:
	/* get mutant space */
	MutantSpace & get_space() const { return mspace; }

	/* whether specified operator refers to some group */
	bool has_operator(const std::string & op) const { return operator_groups.count(op) > 0; }
	/* get the group of specified operator name */
	const MutGroup & get_operator_group(const std::string &) const;
	/* get the set of operator names in mutants */
	const std::set<std::string> & get_operators() const { return operators; }
	/* get group for global subsuming mutants */
	const MutGroup & get_global_group() const;

	/* add, end, clear */
	friend class MutAnalyzer;

protected:
	/* create groups by level in given space */
	MutLevel(MutantSpace & space) : mspace(space), global_group(nullptr), operators(), operator_groups() {}
	/* deconstructor */
	~MutLevel() { clear(); }
	/* add mutant into the level-data */
	void add(Mutant::ID, const BitSeq &);
	/* classify the mutants into specified level */
	void end();
	/* clear all the groups in the level */
	void clear();

protected:
	/* whether the operator is valid name for experiment */
	bool valid_operator(const std::string &);
	/* get the mutant id in cluster */
	Mutant::ID get_mutant_of(const MuCluster &);
private:
	/* mutant space */
	MutantSpace & mspace;
	/* group for global subsumings */
	MutGroup * global_group;
	/* set of operator names */
	std::set<std::string> operators;
	/* map from operator to the group of subsuming mutants */
	std::map<std::string, MutGroup *> operator_groups;
};
/* analyzer for subsuming mutants */
class MutAnalyzer {
public:
	/* create analyzer for subsuming experiment */
	MutAnalyzer(MutantSpace & space) : data(space) {}
	/* deconstructor */
	~MutAnalyzer() {}

	/* update the data */
	void update(ScoreProducer & producer, ScoreConsumer & consumer) {
		ScoreVector * vector;
		
		data.clear();
		while ((vector = producer.produce()) != nullptr) {
			data.add(vector->get_mutant(), vector->get_vector());
			consumer.consume(vector);
		}
		data.end();
	}
	/* get the experiment data */
	MutLevel & get_data() { return data; }

private:
	MutLevel data;
};
/* outputter for subsuming mutations */
class MutOutputer {
public:
	/* create outputter at specified directory */
	MutOutputer(const File &);
	/* deconstructor */
	~MutOutputer() {}

	/* write ../analysis/mutants.txt, ../analysis/summary.txt ../analysis/ */
	void write(const MutLevel &);

protected:
	/* [Mutants; Tests; Killed; Equivaence; Subsuming;] */
	void write_summary(const MutLevel &, std::ostream &);
	/* [mid, oprt, line, origin, replace, category-op, category] */
	void write_mutants(const MutLevel &, std::ostream &);
	/* [op, mutants, equivalence, op-subsuming, glb-subsuming] */
	void write_distribution(const MutLevel &, std::ostream &);

private:
	File * dir;	/* ../analysis/ */
	/* delete lines and \t from the string text */
	void trim(std::string &);
};

/* TestMachine provides following functions:
	1. select subsuming mutants from selective operators;
	2. generate minimal test set for the selected mutants;
	3. evaluate dominator score of tests in the context.
*/
class TestMachine {
public:
	/* create an unintialized machine */
	TestMachine() : context(nullptr) {}
	/* deconstructor */
	~TestMachine() { close(); }

	/* initial the machine by context of subsuming mutants */
	void start(const MutLevel & data) { close(); context = &data; }
	/* select subsuming mutants from given context */
	void select(const std::set<std::string> &, MutantSet &);
	/* generate minimal test set for given mutants */
	void generate(const MutantSet &, TestSet &);
	/* evaluate the dominator score for tests against subsuming mutants in the context  */
	double evaluate(const TestSet &);
	/* close the machine */
	void close() { context = nullptr; }

private:
	const MutLevel * context;

	/* find the first mutant in given set */
	Mutant::ID find_mutant_of(const MutantSet &);
	/* find the first test in score vector that kill some mutant */
	TestCase::ID find_test_in(const BitSeq &);
	/* whether the cluster can be killed by the tests */
	bool kill_cluster(const MuCluster &, const TestSet &);
};

