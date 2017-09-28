#pragma once

/*
	File : suoprt.h
	Purp : to define interface to access subsuming and equivalent mutants
	Arth : Lin Huan
	Date : July 31th, 2017
*/

// include list
#include "mgraph.h"
#include <stack>

// class declarations
class MuClusterSet;
class OpClusterMap;
class OpMutantMap;
class SOperatorSet;

/* equivalent | subsuming | subsumed clusters */
class MuClusterSet {
public:
	/* classify clusters in MSG by their categories */
	MuClusterSet(const MSGraph &);
	/* deconstructor */
	~MuClusterSet() { 
		subsuming_clusters.clear(); 
		subsumed_clusters.clear();
		eq_cluster = nullptr;
	}

	/* get the subsumption graph */
	const MSGraph & get_graph() const { return graph; }

	/* whether there are equivalent mutants */
	bool has_equivalents() const { return eq_cluster != nullptr; }
	/* get the cluster for equivalent mutants */
	MuCluster & get_equivalents() const;
	/* get the clusters for subsuming mutants */
	const std::set<MuCluster *> & get_subsuming() const { return subsuming_clusters; }
	/* get the clusters for subsumed mutants */
	const std::set<MuCluster *> & get_subsumed() const { return subsumed_clusters; }

	/* get the category of cluster: subsumed, equivalent, subsuming */
	char category_of(MuCluster &) const;

	// constant arguments
	const static char Equivalent = 'e';	/* equivalent category */
	const static char Subsuming = 's';	/* subsuming category */
	const static char Subsumed	= 'r';	/* subsumed category */
	const static char NotBelong = 'n';	/* not-belong category */

private:
	/* mutant subsumption graph as source */
	const MSGraph & graph;
	/* cluster to the equivalent mutants */
	MuCluster * eq_cluster;
	/* subsuming mutants clusters */
	std::set<MuCluster *> subsuming_clusters;
	/* clusters for subsumed mutants */
	std::set<MuCluster *> subsumed_clusters;
};
/* map from operators to their clusters */
class OpClusterMap {
public:
	/* create map from operator to clusters (or vice) */
	OpClusterMap(const MSGraph &);
	/* deconstructor */
	~OpClusterMap();

	/* get the subsumption graph where this map is defined */
	const MSGraph & get_graph() const { return graph; }
	/* get the operators set as key */
	const std::set<std::string> & get_operators() const { return operators; }
	/* get the clusters in MSG as key */
	const std::set<MuCluster *> & get_clusters() const { return clusters; }

	/* whether specified operator is defined in map */
	bool has_operator(const std::string & op) const { return operators.count(op) > 0; }
	/* whether the specified cluster is defined */
	bool has_cluster(MuCluster & cluster) const { return clusters.count(&cluster) > 0; }
	/* get the clusters referred by this operators */
	const std::set<MuCluster *> & get_clusters_of(const std::string &) const;
	/* get the operators referred by this cluster */
	const std::set<std::string> & get_operators_of(MuCluster &) const;

private:
	/* subsumption graph as source */
	const MSGraph & graph;
	/* operators key set */
	std::set<std::string> operators;
	/* clusters key set  */
	std::set<MuCluster *> clusters;
	/* map from cluster to operators */
	std::map<MuCluster *, std::set<std::string> *> cop_map;
	/* map from operator to clusters */
	std::map<std::string, std::set<MuCluster *> *> opc_map;
};
/* map from operator to mutants */
class OpMutantMap {
public:
	/* create map from mutant to operator */
	OpMutantMap(MutantSpace &);
	/* deconstructor */
	~OpMutantMap() {
		auto beg = op_map.begin();
		auto end = op_map.end();
		while (beg != end)
			delete (beg++)->second;
		op_map.clear(); operators.clear();
	}

	/* get the mutant space as source */
	MutantSpace & get_space() const { return mspace; }
	/* get the set of operators in mutant space */
	const std::set<std::string> & get_operators() const { return operators; }
	/* whether there is mutant to this operator */
	bool has_operator(const std::string & op) const { return operators.count(op) > 0; }
	/* get the mutants id set for specified operator */
	const std::set<Mutant::ID> & get_mutants_of(const std::string &) const;

private:
	MutantSpace & mspace;
	std::set<std::string> operators;
	std::map<std::string, std::set<Mutant::ID> *> op_map;
};
/* to identify subsuming operators */
class SOperatorSet {
public:
	/* create set for subsuming operators */
	SOperatorSet(const MuClusterSet &, const OpClusterMap &, const OpMutantMap &);
	/* deconstructor */
	~SOperatorSet();

	/* get the clusters set object */
	const MuClusterSet & get_clusters() const { return clusters; }
	/* get the map from operator to clusters */
	const OpClusterMap & get_mappings() const { return mappings; }
	/* get the map from operators to mutants */
	const OpMutantMap & get_operators() const { return operators; }
	/* get the subsuming operators (100-sufficient) */
	const std::set<std::string> & get_subsuming_operators() const { return SOPSet; }
	/* get the alpha-sufficient operators */
	bool get_sufficient_operators(double alpha, std::set<std::string> & ops) {
		if (alpha < 0 || alpha > 1) {
			CError error(CErrorType::InvalidArguments, 
				"SOperatorSet::get_sufficient_operators", 
				"Invalid alpha: " + std::to_string(alpha));
			CErrorConsumer::consume(error);
			exit(CErrorType::InvalidArguments);
		}
		else return solve_by_bounds(alpha, ops);
	}
	/* get the alpha-sufficient operators by greedy algorithm */
	bool get_sufficient_operators_fast(double alpha, std::set<std::string> & ops) {
		if (alpha < 0 || alpha > 1) {
			CError error(CErrorType::InvalidArguments,
				"SOperatorSet::get_sufficient_operators",
				"Invalid alpha: " + std::to_string(alpha));
			CErrorConsumer::consume(error);
			exit(CErrorType::InvalidArguments);
		}
		else return solve_by_greedy(alpha, ops);
	}

private:
	/* set of clusters (by categories) */
	const MuClusterSet & clusters;
	/* map from operators to clusters */
	const OpClusterMap & mappings;
	/* map from operator to mutants */
	const OpMutantMap & operators;
	/* set of subsuming operators set */
	std::set<std::string> SOPSet;

	/* type for state of search item */
	typedef struct {
		std::vector<std::string> op_state;	/* operator for selection */
		std::set<MuCluster *> sc_state;		/* subsuming mutants used */
		int cursor;			/* index to the next operator */
		unsigned value;		/* evaluation hold by current state */
	} _Item;

	std::stack<_Item *> _stack;	/* stack for tree-iteration */

protected:
	/* solve the operators by greedy algorithm */
	bool solve_by_greedy(double alpha, std::set<std::string> &);
	/* solve the operators by bound-limit algorithm */
	bool solve_by_bounds(double alpha, std::set<std::string> &);
	
	/* pop the top item at stack */
	void pope_item();

	/* initialize the state items by initial seed */
	void bound_init_item(const std::set<std::string> &);
	/* push the next item by select the operator at parent and update requirements of covered cluster */
	void bound_push_item();
	/* generate solutions directly from the top item in the stack */
	void top_solutions(std::set<std::string> &);

	/* initialize the stack items by greedy algorithm */
	void greedy_init_item(const std::set<MuCluster *> &);
	/* push the solution for stack of greedy algorithm */
	void greedy_push_item();
	/* generate solutions by the tagged eliminated operator in the lines of stack */
	void lin_solutions(std::set<std::string> &);

	/* get the number of mutants created by this operator */
	unsigned evaluate(const std::string &) const;
};
/* writer for subsuming operators */
class SuOprtWriter {
public:
	/* create writer for subsuming operators */
	SuOprtWriter(const TestSpace & tset) : dir(nullptr), tspace(tset) {}
	/* deconstructor */
	~SuOprtWriter() { close(); }

	/* open ../analysis */
	void open(const File & root) {
		close();	/* initialization */

		/* get the directory */
		const std::vector<File *> & files = root.list_files();
		for (int i = 0; i < files.size(); i++) {
			File * file = files[i];
			if (file->get_local_name() == "analysis") {
				dir = file; break;
			}
		}

		/* validate */
		if (dir == nullptr) {
			CError error(CErrorType::InvalidArguments, 
				"SuOprtWriter::open", 
				"\"analysis\" is not found at: " + root.get_path());
			CErrorConsumer::consume(error);
			exit(CErrorType::InvalidArguments);
		}
	}
	/* writer information in subsuming operators to ../analysis */
	void write(SOperatorSet &);

	/* generate evaluation on specified alpha-SMOs [alpha, oplist, cov, score, costs] */
	void eval_operators(SOperatorSet &, const std::vector<double> &);

	/* close writer */
	void close() { dir = nullptr; }

private:
	const TestSpace & tspace;
	File * dir;	// ../analysis/

protected:
	/* ../analysis/summary.txt */
	void write_summary(SOperatorSet &, std::ostream &);
	/* ../analysis/coverage.txt */
	void write_coverage(SOperatorSet &, std::ostream &);
	/* ../analysis/score_line.txt */
	void write_scoreln(SOperatorSet &, std::ostream &);
	/* ../analysis/mutantset.txt */
	void write_mutants(SOperatorSet &, std::ostream &);
	/* ../analysis/operator_score.txt */
	void write_op_scores(SOperatorSet &, std::ostream &);

	/* {coverage, dom_score} */
	void gen_cov_score(SOperatorSet &, const std::set<std::string> &, std::ostream &, TestSet &);
	/* select operators set based on bit-string */
	void select_operators(const std::vector<std::string> &, const BitSeq &, std::set<std::string> &);

	/* get the covered subsuming clusters by specified operator */
	void get_coverage(SOperatorSet &, std::set<MuCluster *> &, const std::string &);
	/* {alpha, oplist, cop, dscore} */
	void write_oplist(SOperatorSet &, const std::set<std::string> &, std::ostream &, TestSet &);
	/* remove spaces (\t, \n) in text */
	void trim_spaces(std::string &);
};
/* machine for evaluating selective operators */
class TestMachine {
public:
	/* create a closed machine for test generation */
	TestMachine(const SOperatorSet & ctxt) : context(ctxt) {}
	/* deconstructor */
	~TestMachine() {}

	/* generate a minimal test set for subsuming mutants in given operators */
	void generate_by_operators(TestSet &, const std::set<std::string> &);
	/* generate test set from clusters requirement */
	void generate_by_requirement(TestSet &, const std::set<MuCluster *> &);
	/* evaluate the dominator score of given test suite */
	double evaluate(const TestSet &);

	/* get the context where tests are generated */
	const SOperatorSet & get_context() const { return context; }

private:
	/* context for test generation */
	const SOperatorSet & context;

	/* generate a random index in [0, n - 1] */
	BitSeq::size_t gen_random_seed(BitSeq::size_t);
	/* find the kth 1 in bit sequence*/
	TestCase::ID find_test_at(const BitSeq &, TestCase::ID);
	/* whether the tests can kill specified cluster */
	bool is_killed(const TestSet &, const MuCluster &);

protected:
	/* generate tests by greedy algorithm */
	void greedy_generate_tests(TestSet &, const std::set<MuCluster *> &);
};

