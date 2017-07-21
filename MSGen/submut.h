#pragma once

/*
	-File: submut.h
	-Date: July 6th 2017
	-Purp: to define interface for subsuming mutants and analysis
	-Clas:
		[1] class MutGroup
*/

#include "mgraph.h"
#include <stack>

/* class declaration */
class MutGroup;
class MutLevel;
class MutAnalyzer;
class OperatorWriter;

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

/* algorithm to search the minimal subsuming operators for mutants */
class SuOperatorSearch {
public:
	/* create a searcher for identifying subsuming operator */
	SuOperatorSearch(const MutLevel & data) 
		: context(data), op_scs(), sc_ops(), _stack() {}
	/* deconstructor */
	~SuOperatorSearch() { finish(); }

	/* initialize the searcher state */
	void start();
	/* identify the minimal subsuming operator by branch algorithm */
	void solve(std::set<std::string> &);
	/* clear the searcher state machine */
	void finish();

private:

	/* type for state of search item */
	typedef struct {
		std::vector<std::string> op_state;	/* operator for selection */
		std::set<MuCluster *> sc_state;		/* subsuming mutants used */
		int cursor;			/* index to the next operator */
		unsigned value;		/* evaluation hold by current state */
	} _Item;

	const MutLevel & context;		/* context of problem */
	
	std::map<std::string, std::set<MuCluster *> *> op_scs;	/* map from operator to its cover mutants */
	std::map<MuCluster *, std::set<std::string> *> sc_ops;	/* map from mutant to its cover operators */
	std::stack<_Item *> _stack;					/* stack for solution states */

	/* get the cluster equals with source in SC */
	MuCluster * belong_to(MuCluster *, const std::set<MuCluster *> &);
	/* generate solution from stack items */
	void gen_solution(std::set<std::string> &);

	/* evaluate the value of operator */
	unsigned eval(const std::string &);
	/* evaluate operator by the size of generated mutants */
	unsigned eval_mutants(const std::string & op);
	/* evaluate operator by its number  */
	unsigned eval_operats(const std::string & op);
	/* evaluate operator by the size of op-subsuming mutants */
	unsigned eval_smutant(const std::string & op);

protected:
	/* {op_scs.keys; sc_ops.keys; -1; 0;} */
	void int_stat();
	/* push the next item based on previous top item state and selected index */
	void psh_stat();
	/* pop the top item and delete its space */
	void pop_stat();
};
/* writer to identify subsuming operators and compute the contributions */
class OperatorWriter {
public:
	/* create writer for subsuming operators */
	OperatorWriter() : dir(nullptr) {}
	/* deconstructor */
	~OperatorWriter() { close(); }

	/* open ../analysis/ for writing results */
	void open(const File &);
	/* write summary.txt and distribution.txt */
	void write(MutLevel &, const std::set<std::string> &);
	/* close the writer */
	void close();

protected:
	/* {Mut, Equiv, Oprt, SOp, SOp-Mut, SuM} */
	void write_summary(MutLevel &, const std::set<std::string> &, std::ostream &);
	/* {op, Mut, Equiv, SMut, Cop, Eop} */
	void write_contribution(MutLevel &, std::ostream &);

private:
	/* ../analysis/ */
	const File * dir;
	/* whether cluster equals with one of the candidate */
	MuCluster * belong_to(const MuCluster &, const std::set<MuCluster *> &);
};

/* machine for evaluating selective operators */
class TestMachine {
public:
	/* create a closed machine for test generation */
	TestMachine(const MutLevel & ctxt) : context(ctxt), operators() {}
	/* deconstructor */
	~TestMachine() { close(); }

	/* set the context for test generation */
	void start(const std::set<std::string> &);
	/* generate a minimal test set for subsuming mutants in given operators */
	void generate(TestSet &);
	/* evaluate the dominator score of given test suite */
	double evaluate(const TestSet &);
	/* close the test machine for generation */
	void close() { operators.clear(); }

	/* get the context where tests are generated */
	const MutLevel & get_context() const { return context; }
	/* get the subsuming operators selected by start() */
	const std::set<std::string> & get_selected_operators() const { return operators; }

private:
	/* context for test generation */
	const MutLevel & context;
	/* selected operators */
	std::set<std::string> operators;

protected:
	/* select minimal test template for given requirements (greedy algorithm), the requirements will be */
	void select_minimal_template(
		const std::set<MuCluster *> & requirements,
		std::vector<BitSeq *> & templates,
		const TestSpace & test_space
	);
	/* find the kth 1 in bit sequence*/
	TestCase::ID find_test_at(const BitSeq &, TestCase::ID);
	/* generate tests from template by given seeds */
	void generate_test_suite(
		const std::vector<BitSeq *> & templates,
		const std::vector<unsigned> & seeds,
		TestSet & tests
	);
	/* whether the tests can kill specified cluster */
	bool is_killed(const TestSet &, const MuCluster &);
};