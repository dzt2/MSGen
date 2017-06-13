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
	TypedMutantSet(MSGraph &);
	~TypedMutantSet();

	const MutantSet & get_mutants() const;
	const MutantSet & get_stubborn_mutants() const;
	const MutantSet & get_subsuming_mutants() const;
	const MutantSet & get_subsumed_mutants() const;



};










