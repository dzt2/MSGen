#pragma once

/*
	-File : submut.h
	-Arth : Lin Huan
	-Date : June 12th, 2017
	-Purp : to define interface to compute and analyze subsuming mutations 
	-Clas :
		[1] TypedMutantSet
		[2] 
*/

#include "mgraph.h"

class TypedMutantSet;

/* Mutants in MSG are classified in three types:
	I.	stubborn mutants: cannot be killed by any tests in finite T;
	II. subsuming mutants: directly subsumed by stubborn mutants, i.e. non-redundant set of mutations;
	III. subsumed mutants: subsumed by subsuming mutants, i.e. redundant mutants.
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
	/* get the subsuming mutants (non-redundant) */
	const MutantSet & get_subsuming_mutants() const { return *subsuming_mutants; }
	/* get the subsumed mutants (redundant) */
	const MutantSet & get_subsumed_mutants() const { return *subsumed_mutants; }

	/* get the mutant subsumption graph that defines the set*/
	MSGraph & get_graph() const { return graph; }
	/* get the cluster of stubborn mutants */
	MuCluster & get_stubborn_cluster() const { return *stubborn_cluster; }
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
	/* subsuming mutants (non-redundant) */
	MutantSet * subsuming_mutants;
	/* subsumed mutants (redundant) */
	MutantSet * subsumed_mutants;

	/* stubborn cluster (unkilled) */
	MuCluster * stubborn_cluster;
	/* subsuming clusters (minimal non-redundant) */
	std::set<MuCluster *> subsuming_clusters;
	/* subsumed clusters (redundant) */
	std::set<MuCluster *> subsumed_clusters;
};
