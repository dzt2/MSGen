#pragma once

/*
	-File : submut.h
	-Arth : Lin Huan
	-Date : June 12th, 2017
	-Purp : to define interface to compute and analyze subsuming mutations 
	-Clas :
		[1] TypedMutantSet
		[2] TypedOutputter
		[4] GreedyTester
*/

#include "mgraph.h"

class TypedMutantSet;
class TypedOutputter;
class GreedyTester;

/* Mutants in MSG are classified in three types:
	I.	stubborn mutants: cannot be killed by any tests in finite T;
	II. subsuming mutants: directly subsumed by stubborn mutants, i.e. non-redundant set of mutations;
	III. subsumed mutants: subsumed by subsuming mutants, i.e. redundant mutants.
	IV. H0-hard mutants: mutants in H[1] (with least non-zero degree)
*/
class TypedMutantSet {
public:
	/* create a typed-set of mutations by MSG */
	TypedMutantSet(MSGraph &);
	/* deconstructor */
	~TypedMutantSet();

	/* get all the mutants in this set */
	const MutantSet & get_mutants() const { return mutants; }
	/* get the stubborn mutants */
	const MutantSet & get_stubborn_mutants() const { return *stubborn_mutants; }
	/* get the hard mutants (killed by least tests) */
	const MutantSet & get_hard_mutants() const { return *hard_mutants; }
	/* get the subsuming mutants (non-redundant) */
	const MutantSet & get_subsuming_mutants() const { return *subsuming_mutants; }
	/* get the subsumed mutants (redundant) */
	const MutantSet & get_subsumed_mutants() const { return *subsumed_mutants; }

	/* get the mutant subsumption graph that defines the set*/
	MSGraph & get_graph() const { return graph; }
	/* get the cluster of stubborn mutants */
	MuCluster & get_stubborn_cluster() const { return *stubborn_cluster; }
	/* get the clusters for hard mutants (killed by least tests) */
	const std::set<MuCluster *> & get_hard_clusters() const { return hard_clusters; }
	/* get the clusters of subsuming mutants */
	const std::set<MuCluster *> & get_subsuming_clusters() const { return subsuming_clusters; }
	/* get the clusters for subsumed mutants */
	const std::set<MuCluster *> & get_subsumed_clusters() const { return subsumed_clusters; }

private:
	/* subsumption graph */
	MSGraph & graph;
	/* set of all mutations */
	const MutantSet & mutants;

	/* stubborn mutants (unkilled) */
	MutantSet * stubborn_mutants;
	/* hard mutants are killed by least tests */
	MutantSet * hard_mutants;
	/* subsuming mutants (non-redundant) */
	MutantSet * subsuming_mutants;
	/* subsumed mutants (redundant) */
	MutantSet * subsumed_mutants;

	/* stubborn cluster (unkilled) */
	MuCluster * stubborn_cluster;
	/* set of clusters for hard mutants */
	std::set<MuCluster *> hard_clusters;
	/* subsuming clusters (minimal non-redundant) */
	std::set<MuCluster *> subsuming_clusters;
	/* subsumed clusters (redundant) */
	std::set<MuCluster *> subsumed_clusters;
};
/* to output the statistical information for typed mutants */
class TypedOutputter {
public:
	/* create outputter for typed mutants */
	TypedOutputter() : dir(nullptr) {}
	/* deconstructor */
	~TypedOutputter() { close(); }

	/* set the ../ where defines director of "analysis" for storing generated files */
	void open(const File &);
	/* generate ../analysis/class.txt from typed mutants */
	void output_classification(const TypedMutantSet &);
	/* generate ../analysis/summary.txt */
	void output_summary(const TypedMutantSet &);
	/* output ../analysis/graph.txt */
	void output_graph(const MSGraph &);

	/* close the outputter */
	void close() { dir = nullptr; }

private:
	/* ../analysis */
	File * dir;
	/* remove all '\n' from the code */
	void trim_lines(std::string &);
	/* output category for mutant by operator name */
	void output_categories(const std::string &, std::string &, std::string &);

protected:
};
