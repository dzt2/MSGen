#pragma once

/*
	-File : smerge.h
	-Arth : Lin Huan
	-Date : June 26th 2017
	-Purp : to define interfaces for deriving subsuming mutants step by step.
	-Clas : 
		[1] SuMutantSet
		[2] SuMutantMerger
*/

#include "mgraph.h"

class SuMutantSet;
class SuMutantMerger;

/* interface to access subsuming mutants from MSG */
class SuMutantSet {
public:
	/* create set for subsuming mutants */
	SuMutantSet(MSGraph & g) : graph(g) { 
		mutants = g.get_space().create_set(); 
		update_mutants();
	}
	/* deconstructor */
	~SuMutantSet() { graph.get_space().delete_set(mutants); }

	/* get the graph that defines the subsuming mutants */
	MSGraph & get_graph() const { return graph; }
	/* get the set of subsuming mutants in graph */
	const MutantSet & get_subsuming_mutants() const { return *mutants; }

	/* whether the mutant is in the set */
	bool has_mutant(Mutant::ID mid) const { return mutants->has_mutant(mid); }
	/* get cluster of specified mutant */
	MuCluster & get_cluster_of(Mutant::ID mid) const { return graph.get_cluster_of(mid); }

	/* derive subsuming mutants from graph */
	void update_mutants();

private:
	MSGraph & graph;
	MutantSet * mutants;
};
/* to merge the subsuming mutants to gain new subsuming set */
class SuMutantMerger {
public:
	/* merger for subsuming mutants between several graphs */
	SuMutantMerger() : builder() {}
	/* deconstructor */
	~SuMutantMerger() { close(); }

	/* the subsuming mutants (between graphs) are maintained in this set */
	void open(SuMutantSet &);
	/* add the set of subsuming mutants for filtering */
	void append(SuMutantSet &);
	/* extract the subsuming mutants from child graphs */
	void extract();
	/* close the merger for next computation */
	void close() { builder.close(); answer = nullptr; }

private:
	MSGBuilder builder;
	SuMutantSet * answer;
};

/* core data space for experiments of subsuming mutants */
class SuMutantExperimentCore {
public:
	/* create the core data for experiment */
	SuMutantExperimentCore(MutantSpace & space) 
		: mspace(space), keys(), mut_map() {}
	/* deconstructor */
	~SuMutantExperimentCore() { clear_mutants(); }

	/* get the mutant space */
	MutantSpace & get_space() const { return mspace; }

	/* get the keys for accessing subsuming mutants */
	const std::set<std::string> & get_keys() const { return keys; }
	/* whether there are subsuming mutants for specified operator  */
	bool has_subsuming_mutants(const std::string & key) const { return (mut_map.count(key) > 0); }
	/* get the subsuming mutants for specified operator */
	const MutantSet & get_subsuming_mutants(const std::string &) const;

	/* get_mutants | clear_mutants */
	friend class SuMutantExperimentDriver;
protected:
	/* get the subsuming mutants for specified key 
		1) if not defined, create a new object;
		2) otherwise, return the existing one.
	*/
	SuMutantSet & get_mutants(const std::string &);
	/* clear the mutants in the data core */
	void clear_mutants();

private:

	MutantSpace & mspace;
	/* keys for operator names */
	std::set<std::string> keys;
	/* map from names to subsuming mutations */
	std::map<std::string, SuMutantSet *> mut_map;
};
/* experiment driver */
class SuMutantExperimentDriver {
public:
	/* driver for experiment */
	SuMutantExperimentDriver() : core(nullptr) {}
	/* deconstructor */
	~SuMutantExperimentDriver() { finish(); }

	/* start up the platform of experiment */
	void start(SuMutantExperimentCore & c) { finish(); core = &c; }
	/* derive the subsuming mutants by their operator, directly */
	void derive_operator_I(ScoreProducer &, ScoreConsumer &);
	/* derive the subsuming mutants by grouped operators (designed) */
	void derive_operator_II();
	/* derive the global subsuming mutants for all the mutants */
	void derive_global_III();
	/* terminate the experiment */
	void finish() { core = nullptr; }

private:
	SuMutantExperimentCore * core;
};