#pragma once

/*
-File : ctrace.h
-Arth : Lin Huan
-Date : May 15th 2017
-Purp : to define interfaces for coverage
-Clas :
[1] LineCoverage
[2] FileCoverage
[3] CoverageSpace
[4] CoverageSource
[5] CoverageLoader
[6] CTrace
*/

#include "bitseq.h"
#include "cerror.h"
#include "cfile.h"
#include "ctest.h"
#include <map>

// declarations 
class LineCoverage;
class FileCoverage;
class CoverageSpace;
class CoverageSource;
class CoverageLoader;
class CTrace;

/* to represent tests that cover the specific line in code file */
class LineCoverage {
protected:
	/* construct a line coverage in code file */
	LineCoverage(const FileCoverage & fil, size_t lin, BitSeq::size_t tnum)
		: file(fil), line(lin), vec(tnum) {}
	/* deconstructor */
	~LineCoverage() {}

	/* set the kth bit in bit-string as 1 */
	void cover(const BitSeq::size_t);
public:
	/* get the file coverage of this line */
	const FileCoverage & get_file() const { return file; }
	/* get the line of coverage in code file, starting from 1*/
	const size_t & get_line() const { return line; }
	/* get the bit-string that represents  */
	const BitSeq & get_vector() const { return vec; }

	/* create and delete and set coverage */
	friend class FileCoverage;
private:
	/* file coverage where this line is defined */
	const FileCoverage & file;
	/* number of line to be covered, starting from 1 */
	const size_t line;
	/* bitstring to represent the tests that covers this line */
	BitSeq vec;
};
/* coverage of a code file */
class FileCoverage {
protected:
	/* constructor */
	FileCoverage(const CoverageSpace &, const CodeFile &, const TestSet &);
	/* deconstructor */
	~FileCoverage() { release(); }

	/* load the lines in code file, and update the map from bit-string index to test case id at the same time! */
	void load_in();
	/* clear the lines, delete all its line coverage in it, and clear test-template at the same time */
	void release();
	/* set kth line is covered by specified test case */
	void cover(size_t, TestCase::ID);

public:
	/* get the space where the file coverage is defined */
	const CoverageSpace & get_space() const { return space; }
	/* get the code file of this coverage */
	const CodeFile & get_file() const { return cfile; }
	/* get the set of tests of this file coverage */
	const TestSet & get_tests() const { return tests; }

	/* whether the kth line is covered */
	bool is_covered(const size_t) const;
	/* get the coverage of kth line */
	const LineCoverage & get_line_coverage(const size_t) const;

	/* create, delete */
	friend class CoverageSpace;
	/* load-release and set-cover */
	friend class CoverageLoader;
private:
	/* space where this file coverage is defined */
	const CoverageSpace & space;
	/* code file of this coverage */
	const CodeFile & cfile;
	/* set of tests to cover this file */
	const TestSet & tests;
	/* list that refers to coverage in each line */
	std::vector<LineCoverage *> lines;
};
/* space for all coverage information */
class CoverageSpace {
protected:
	/* construct a coverage space of all code files in code space */
	CoverageSpace(const CTrace &, const CoverageSource &, const CodeSpace &, const TestSpace &);
	/* deconstructor */
	~CoverageSpace() { this->clear(); }
public:
	/* get the project where this space is defined */
	const CTrace & get_project() const { return project; }
	/* get the source of the coverage space */
	const CoverageSource & get_source() const { return source; }
	/* get the code space for reference */
	const CodeSpace & get_code_space() const { return code_space; }
	/* get the test space for reference */
	const TestSpace & get_test_space() const { return test_space; }

	/* whether there has file coverage for code file */
	bool has_file_coverage(const CodeFile &) const;
	/* get the coverage of specified code file */
	FileCoverage & get_file_coverage(const CodeFile &) const;
	/* set the specified code file as covered by specified test set, this will clear original file coverage! */
	void add_file_coverage(const CodeFile &, const TestSet &);
	/* remove the coverage of specified code file from space */
	void del_file_coverage(const CodeFile &);
	/* clear all the coverage in the space */
	void clear();

	// create and delete
	friend class CTrace;
private:
	/* project of coverage */
	const CTrace & project;
	/* coverage source */
	const CoverageSource & source;
	/* code space for reference */
	const CodeSpace & code_space;
	/* test space for reference */
	const TestSpace & test_space;
	/* map from code file to their file coverage */
	std::map<const CodeFile *, FileCoverage *> files;
};
/* source of ../traces/ */
class CoverageSource {
protected:
	/* constructor */
	CoverageSource(const File & d) : dir(d) {}
	/* deconstructor */
	~CoverageSource() {}
public:
	/* get ../traces/ */
	const File & get_trace_dir() const { return dir; }

	friend class CTrace;
private:
	/* ../traces/ */
	const File & dir;
};
/* to load the lines in code file and extract coverage from ../trace */
class CoverageLoader {
protected:
	CoverageLoader() {}
	~CoverageLoader() {}

	/* load the lines in code file to the file coverage */
	void load_in(FileCoverage & fcov) const { fcov.load_in(); }
	/* release all lines in file coverage */
	void release(FileCoverage & fcov) const { fcov.release(); }
	/* parse the xxx.c.gcov to coverage information in file coverage */
	void parse(FileCoverage &) const;

public:
	// create and delete and load 
	friend class CTrace;

private:
	/* set the coverage of test in specified txx for .gcov file */
	void parse(FileCoverage &, TestCase::ID, const File &) const;
};
/* project for trace */
class CTrace {
public:
	/* construct project for trace analysis */
	CTrace(const File &, const CodeSpace &, const TestSpace &);
	/* deconstructor */
	~CTrace() { delete space; delete source; }

	/* get ../ */
	const File & get_root() const { return dir; }
	/* get the coverage space of this project */
	CoverageSpace & get_space() const { return *space; }
	/* load the coverage in file-coverage by ../traces/xxx.c/ */
	bool load_coverage(const CodeFile &);

private:
	const File & dir;
	CoverageSpace * space;
	CoverageSource * source;
	CoverageLoader loader;
};
