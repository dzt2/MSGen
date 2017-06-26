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

	/* type for category */
	typedef char Category;
	/* category for equivalent mutant in this set */
	static const char equivalent = 0;
	/* category for subsuming mutant in this set */
	static const char subsuming = 1;
	/* category for subsumed mutant in this set */
	static const char subsumed = 2;
	/* category that describe mutants out of this set */
	static const char not_belong = 3;

	/* create set for subsuming mutants */
	SuMutantSet(MSGraph & g) : graph(g) { 
		mutants = g.get_space().create_set(); 
		update_mutants();
	}
	/* deconstructor */
	~SuMutantSet() { graph.get_space().delete_set(mutants); }

	/* get the graph that defines the subsuming mutants */
	MSGraph & get_graph() const { return graph; }

	/* get the equivalent cluster */
	MuCluster * get_equivalent_cluster() const { return equivalent_cluster; }
	/* get the set of clusters for subsuming mutants */
	const std::set<MuCluster *> & get_subsuming_clusters() const { return clusters; }

	/* get the set of equivalent mutants, if there are no equivalent mutants in the set, return null */
	const MutantSet * get_equivalent_mutants() const {
		if (equivalent_cluster == nullptr) 
			return nullptr;
		else 
			return &(equivalent_cluster->get_mutants());
	}
	/* get the set of subsuming mutants in graph */
	const MutantSet & get_subsuming_mutants() const { return *mutants; }
	/* get all the mutants in the set */
	const MutantSet & get_all_mutants() const { return graph.get_mutants(); }

	/* whether the mutant is in the set */
	bool has_mutant(Mutant::ID mid) const { return mutants->has_mutant(mid); }
	/* get cluster of specified mutant */
	MuCluster & get_cluster_of(Mutant::ID mid) const { return graph.get_cluster_of(mid); }

	/* get the category of the specified mutant */
	Category get_category_of(Mutant::ID mid) const {
		if (graph.get_mutants().has_mutant(mid)) {
			if (equivalent_cluster != nullptr
				&& equivalent_cluster->has_mutant(mid))
				return equivalent;
			else if (mutants->has_mutant(mid))
				return subsuming;
			else return subsumed;
		}
		else not_belong;
	}

	/* derive subsuming mutants from graph */
	void update_mutants();

private:
	/* graph that defines the set of subsuming mutants */
	MSGraph & graph;

	/* set of subsuming mutants */
	MutantSet * mutants;

	/* cluster for equivalent mutants */
	MuCluster * equivalent_cluster;
	/* clusters for subsuming mutants */
	std::set<MuCluster *> clusters;
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
	const SuMutantSet & get_subsuming_mutants(const std::string &) const;

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
	/* core data for experiment platform */
	SuMutantExperimentCore * core;
};
/* outputter for experiment result */
class SuMutantExperimentOutput {
public:
	/* output for experiment results */
	SuMutantExperimentOutput() : dir(nullptr) {}
	/* deconstructor */
	~SuMutantExperimentOutput() { close(); }

	/* output files in this directory */
	void open(const File & d) { close(); dir = &d; }
	/* output the information of result to target directory */
	void write(const SuMutantExperimentCore &);
	/* close the output stream */
	void close() { dir = nullptr; }

protected:
	/* for _summary_.txt */
	void write_summary(const SuMutantExperimentCore &, std::ostream &);
	/* for operator.txt */
	void write_mutants(const SuMutantSet &, std::ostream &);

private:
	/* directory */ 
	const File * dir;
	/* remove all the spaces in the string */
	void trim(std::string & str);
};