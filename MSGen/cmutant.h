#pragma once

/*
-File : cmutant.h
-Arth : Lin Huan
-Date : May 16th 2017
-Purp : to define interfaces to access mutants of code file
-Clas :
[1] class Mutation		(in one location)
[2] class Mutant		(to corresponding faulty version)
[3] class MutantSpace	(where mutant is defined)
[4] class MutantSource	(../muta/xxx/ for a source/preprocessed file)
[5] class MutantSet		(subset of mutants in space)
[6] class CMutant		(project to manage mutants in ../muta/{*})
*/

#include "bitseq.h"
#include "cfile.h"
#include <map>

// declarations
class Mutation;
class Mutant;
class MutantSpace;
class MutantSet;
class MutantSource;
class CMutant;

/* mutation is a syntactic change at one location in code file */
class Mutation {
protected:
	/* construct a mutation by specifying its location and text to replace original code */
	Mutation(const Mutant & mut, const CodeFile & cfile, size_t bias, size_t len,
		const std::string & rep) : mutant(mut), loc(cfile, bias, len), replace(rep) {}
	~Mutation() {}

public:
	/* get the mutant of this mutation */
	const Mutant & get_mutant() const { return mutant; }
	/* get the location where this mutation is applied */
	const CodeLocation & get_location() const { return loc; }
	/* get the text to replace original code at location */
	const std::string & get_replacement() const { return replace; }

	/* create and delete */
	friend class Mutant;
private:
	/* mutant to this mutation belongs */
	const Mutant & mutant;
	/* location where the mutation is applied */
	const CodeLocation loc;
	/* text to replace original code at location */
	std::string replace;
};
/* mutant consists of several mutations in one code file (higher-order) */
class Mutant {
public:
	typedef unsigned int ID;

	/* get the space where mutant is defined*/
	const MutantSpace & get_space() const { return space; }
	/* get the id of this mutant */
	ID get_id() const { return mid; }
	/* get the name of mutant operator */
	const std::string & get_operator() const { return oprt; }

	/* get the order of mutant */
	size_t get_orders() const { return orders; }
	/* get the kth mutation in this mutant */
	const Mutation & get_mutation(size_t) const;

	/* create and delete */
	friend class MutantSpace;
	/* set mutation */
	friend class MutantLoader;

protected:
	/* construct a mutant in space with specified degree */
	Mutant(const MutantSpace &, ID, const std::string &);
	/* deconstructor */
	~Mutant() { clear_mutations(); }

	/* set the order of mutant, i.e. number of mutations, this will clear original mutations from list */
	void set_orders(size_t);
	/* set the mutation at specified location of code file with replacement */
	void set_mutation(size_t, size_t, size_t, const std::string &);
	/* remove all mutations from the mutant */
	void clear_mutations();
private:
	/* space where the mutant is defined */
	const MutantSpace & space;
	/* id of this mutant */
	ID mid;
	/* name of mutant operator */
	std::string oprt;
	/* orders of the mutant, usually 1--2 for Proteum */
	size_t orders;
	/* list of mutations of this mutant */
	Mutation ** mutations;
};
/* source of mutant, ../muta/xxx */
class MutantSource {
protected:
	/* create source of ../muta/xxx/ */
	MutantSource(const File &);
	/* deconstructor */
	~MutantSource() {}

public:
	/* get ../muta/xxx/ */
	const File & get_dir() const { return dir; }
	/* get ../muta/xxx/__xxx.c */
	const File & get_code() const { return *code; }
	/* get ../muta/xxx/mschema.txt */
	const File & get_schema() const { return *schema; }

	/* create and delete */
	friend class MutantSpace;
private:
	/* ../muta/xxx/ */
	const File & dir;
	/* ../muta/xxx/__xxx.c */
	const File * code;
	/* ../muta/xxx/mschema.txt */
	const File * schema;
};
/* space for ../muta/xxx/ referring to one code file */
class MutantSpace {
protected:
	/* construct a mutant space for specific code file */
	MutantSpace(const CMutant &, const File & dir, const CodeFile &);
	/* deconstructor */
	~MutantSpace() { clear_sets(); clear_mutants(); }

	/* create a new mutant by specifying its operator */
	Mutant * create_mutant(const std::string &);
	/* clear all mutants in the space */
	void clear_mutants();
	/* clear all the mutant sets in the space */
	void clear_sets();
public:
	/* get the project of mutants in space */
	const CMutant & get_project() const { return project; }
	/* get the source of mutant sapce */
	const MutantSource & get_source() const { return source; }
	/* get the code file to which mutants in the space refer */
	const CodeFile & get_code_file() const { return code_file; }

	/* get the number of mutants in space */
	size_t number_of_mutants() const { return mutants.size(); }
	/* get the kth mutant, starting from 0 */
	Mutant & get_mutant(Mutant::ID) const;

	/* get a new mutant set */
	MutantSet * create_set();
	/* whether the mutant set is available in the space */
	bool is_accessible(MutantSet * mset) const { return mutsets.count(mset) > 0; }
	/* delete the mutant set */
	bool delete_set(MutantSet *);

	/* create and delete and load */
	friend class CMutant;
	/* create mutant */
	friend class MutantLoader;
private:
	/* project for mutations */
	const CMutant & project;
	/* ../muta/xxx */
	const MutantSource source;
	/* code file in which mutants are seeded in this space */
	const CodeFile & code_file;
	/* list of mutants */
	std::vector<Mutant *> mutants;
	/* set of allocated set */
	std::set<MutantSet *> mutsets;
};
/* sub-set of mutants in space */
class MutantSet {
protected:
	/* construct an empty set for mutants in the space */
	MutantSet(const MutantSpace & spac) : space(spac),
		vec(space.number_of_mutants()), number(0) {}
	/* deconstructor */
	~MutantSet() {}

public:
	/* get the space of mutants */
	const MutantSpace & get_space() const { return space; }
	/* get the bit-string of mutant set */
	const BitSeq & get_set_vector() const { return vec; }

	/* whether mutant of specified id in space belongs to this set */
	bool has_mutant(Mutant::ID) const;
	/* add a new mutant of space into this set */
	bool add_mutant(Mutant::ID);
	/* delete an existing mutant from this set */
	bool del_mutant(Mutant::ID);
	/* parse the set to its complement */
	bool complement();

	/* get the number of mutants in this set */
	size_t number_of_mutants() const { return number; }

	/* create and delete */
	friend class MutantSpace;
private:
	/* space where mutants are defined */
	const MutantSpace & space;
	/* bit-string to represent */
	BitSeq vec;
	/* number of mutants in set */
	size_t number;
};

/* loading mutant + mutations from mschema.txt for MutantSpace */
class MutantLoader {
protected:
	MutantLoader() {}
	~MutantLoader() {}
	/* load "../muta/xxx/mschema.txt" to the mutant in space
	1) a shallow copy only loads operator and orders of mutant;
	2) a deeping copy also loads mutations into the mutant;
	*/
	void load(MutantSpace &mspace, bool) const;

	friend class CMutant;
private:
	bool get_id(LineReader &, Mutant::ID &) const;
	bool get_order(LineReader &, size_t &) const;
	bool get_operator(LineReader &, std::string &) const;
	bool get_mutations(LineReader &, size_t &, size_t &, std::string &) const;
};
/* project of mutants */
class CMutant {
public:
	/* create mutant project for ../muta/ */
	CMutant(const File &, const CodeSpace &);
	/* deconstructor */
	~CMutant();

	/* get ../muta/ */
	const File & get_muta_dir() const { return *dir; }
	/* get the code space of the mutants */
	const CodeSpace & get_code_space() const { return cspace; }

	/* whether there are mutants for specified code file */
	bool has_mutants(const CodeFile & cfile) const { return spaces.count(&cfile) > 0; }
	/* get the mutant space of specified code file */
	MutantSpace & get_mutants_of(const CodeFile & cfile) const {
		if (spaces.count(&cfile) == 0) {
			CError error(CErrorType::InvalidArguments, "CMutant::get_mutants_of", "Undefined code-file: " + cfile.get_file().get_path());
			CErrorConsumer::consume(error);
		}
		auto iter = spaces.find(&cfile);
		return *(iter->second);
	}
	/* load mutants for space from ../muta/xxx/mschema.txt */
	bool load_mutants_for(MutantSpace & space, bool deep) const {
		const CMutant & project = space.get_project();
		if (&project != this) {
			CError error(CErrorType::Runtime, "CMutant::load_mutants_for", "Undefined space...");
			CErrorConsumer::consume(error);
		}
		loader.load(space, deep); return true;
	}
	/* get the map from code file to mutant space */
	const std::map<const CodeFile *, MutantSpace *> & get_spaces() const { return spaces; }

private:
	const File * dir;
	const CodeSpace & cspace;
	const MutantLoader loader;
	std::map<const CodeFile *, MutantSpace *> spaces;
};

