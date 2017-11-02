#pragma once

/*
-File : ctest.h
-Arth : Lin Huan
-Date : May 15th, 2017
-Purp : to define interfaces for tests.
-Clas :
[1] TestCase
[2] TestSource
[3] TestSpace
[4] TestSet
[5] SuiteType
[6] TestLoader
[7] CTest
*/

#include "bitseq.h"
#include "cfile.h"
#include <map>

class TestCase;
class TestSpace;
class TestSource;
class TestSet;

class TestLoader;
enum TestType;

class CTest;

/* type for program of test suite*/
enum TestType {
	general,
	print_tokens,
	print_tokens2,
	replace,
	schedule,
	schedule2,
	tcas,
	tot_info,

	flex,
	gzip,
	sed,
	space,
};

/* TestCase = {input; output; program;} */
class TestCase {
public:
	/* test ID */
	typedef unsigned int ID;

	/* get the space where test case is defined */
	const TestSpace & get_space() const { return space; }
	/* get the id of this test case in current space */
	const ID & get_id() const { return tid; }
	/* get the input string of this test case */
	const std::string & get_input_string() const { return input_string; }
	/* get the output file (local path) of this test case, in directory of output/ */
	const std::string & get_output_file() const { return output_file; }
	/* get the executionable program to be tested by this test case */
	const ExeFile & get_program() const { return program; }

	friend class TestSpace;
private:
	/* id of test case in current space */
	ID tid;
	/* test space where the test case is defined */
	const TestSpace & space;
	/* program to be executed and tested */
	const ExeFile & program;
	/* string to represent the input content */
	std::string input_string;
	/* local name of test output file, in the context of output directory */
	std::string output_file;
protected:
	/* create a new test case in current space */
	TestCase(const TestSpace & spac, ID id, const ExeFile & prog, const std::string & input)
		: space(spac), tid(id), program(prog), input_string(input) {
		output_file = "t" + std::to_string(tid);
	}
	/* deconstructor */
	~TestCase() {}
};
/* space to define and manage test cases */
class TestSpace {
protected:
	/* create test space on given source */
	TestSpace(const CTest &, const TestSource &, const ExeSpace &);
	/* deconstructor */
	~TestSpace() { clear(); }

public:
	/* get the project of this test space */
	const CTest & get_project() const { return project; }
	/* get the test source */
	const TestSource & get_source() const { return source; }

	/* get the number of tests in the space */
	std::vector<TestCase *>::size_type number_of_tests() const { return test_list.size(); }
	/* get the test case at specified location of list */
	const TestCase & get_test(const TestCase::ID) const;

	/* clear all test cases in the space */
	void clear();
	/* create a new test case in the tail of list (automatically encoded) */
	bool create_test(const std::string &);

	/* get the space of executionable programs */
	const ExeSpace & get_exec_space() const { return exec_space; }

	// create and delete 
	friend class CTest;
private:
	/* project for test space */
	const CTest & project;
	/* source where test data is maintianed */
	const TestSource & source;

	/* list of test cases in the space */
	std::vector<TestCase *> test_list;
	/* set of test inputs in the space */
	std::set<std::string> input_set;

	/* space that refers to directory of "exec/" */
	const ExeSpace & exec_space;
};
/* source of tests to provide reference to data resources */
class TestSource {
protected:
	/* construct test source in given directory */
	TestSource(const File &);
	/* deconstructor */
	~TestSource() {}

public:
	/* get the "inputs/" */
	const File & get_input_dir() const { return *input_dir; }
	/* get the "outputs/" */
	const File & get_output_dir() const { return *output_dir; }
	/* get the "new_outputs/" */
	const File & get_new_output_dir() const { return *new_output_dir; }
	/* get the "suites/" */
	const File & get_suite_dir() const { return *suite_dir; }

	// create and delete 
	friend class CTest;
private:
	/* directory of test inputs */
	const File * input_dir;
	/* directory of test outputs */
	const File * output_dir;
	/* directory of original outputs */
	const File * new_output_dir;
	/* directory of test suite */
	const File * suite_dir;
};
/* to represent a sub-set of all tests in the space by bit-string */
class TestSet {
protected:
	/* create an empty test set in space */
	TestSet(const TestSpace & spac)
		: space(spac), set_vector(spac.number_of_tests()), number(0) {}
	/* deconstructor */
	~TestSet() {}

public:
	/* get the test space for this test set */
	const TestSpace & get_space() const { return space; }
	/* get the bit vector to represent the test set in the space */
	const BitSeq & get_set_vector() const { return set_vector; }

	/* whether the test case of specified id defined in this set */
	bool has_test(TestCase::ID) const;
	/* get the test case in this  */
	const TestCase & get_test(TestCase::ID) const;

	/* clear all the tests in the space */
	void clear();
	/* add a new test into the set */
	bool add_test(TestCase::ID);
	/* delete an existing test from the set */
	bool del_test(TestCase::ID);
	/* parse the set to its complement */
	void complement();
	/* get the number of tests */
	size_t size() const { return number; }

	/* get the set vector */
	const BitSeq & get_vector() const { return set_vector; }

	/* create and delete */
	friend class CTest;

private:
	/* test space in which tests in set are defined */
	const TestSpace & space;
	/* bit string to represent the tests of this set in the space */
	BitSeq set_vector;
	/* number of tests in the set */
	size_t number;
};

/* to parse the test suite files into test cases in space */
class TestLoader {
public:
	/* create a test suite parser */
	TestLoader(const TestSpace &);
	/* deconstructor */
	virtual ~TestLoader();
	/* whether there has next line */
	bool has_next() { return lineptr != nullptr; }
	/* get the next input string (dynamically allocated), if no more inputs, return null */
	virtual std::string * next_inputs() {
		std::string * line = lineptr;
		roll_next(); return line;
	}
protected:
	/* test space */
	const TestSpace & space;

	/* pointer to the next line */
	std::string * lineptr;
	/* reader for file in cur */
	LineReader * reader;
	/* currently accessed file in suites/ */
	std::vector<File *>::const_iterator cur;
	/* end of files in suites/ */
	std::vector<File *>::const_iterator end;

	/* get the next line in suite file(s) */
	void roll_next();
};
/* loader for test suites in print_tokens(2) */
class TestLoader_Printtokens : public TestLoader {
public:
	TestLoader_Printtokens(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Printtokens() {}
	std::string * next_inputs();
};
/* loader for test suites in replace */
class TestLoader_Replace : public TestLoader {
public:
	TestLoader_Replace(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Replace() {}
	std::string * next_inputs();
};
/* loader for test suites in schedule(2) */
class TestLoader_Schedule : public TestLoader {
public:
	TestLoader_Schedule(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Schedule() {}
	std::string * next_inputs();
};
/* loader for test suites in tcas */
class TestLoader_Tcas : public TestLoader {
public:
	TestLoader_Tcas(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Tcas() {}
	std::string * next_inputs();
};
/* parser for test suites in tcas */
class TestLoader_Totinfo : public TestLoader {
public:
	TestLoader_Totinfo(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Totinfo() {}
	std::string * next_inputs();
};
/* parser for flex test suite */
class TestLoader_Flex : public TestLoader {
public:
	TestLoader_Flex(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Flex() {}
	std::string * next_inputs();

private:
	char next_tag(const std::string &, int &);
	bool extracts(const std::string &, int &, std::string &);
	bool extrange(std::string & inlist, std::string &oulist);
};
/* parser for gzip test suite */
class TestLoader_Gzip : public TestLoader {
public:
	TestLoader_Gzip(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Gzip() {}
	std::string * next_inputs();
private:
	char get_tag(const std::string &, int &);
	bool get_content(const std::string &, int &, std::string &);
	bool is_test(const std::string &);
};
/* parser for sed test suite */
class TestLoader_Sed : public TestLoader {
public:
	TestLoader_Sed(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Sed() {}
	std::string * next_inputs();

private:
	char get_tag(const std::string &, int &);
	bool get_content(const std::string &, int &, std::string &);
};
/* parser for space test suite */
class TestLoader_Space : public TestLoader {
public:
	TestLoader_Space(const TestSpace & spac) : TestLoader(spac) {}
	~TestLoader_Space() {}
	std::string * next_inputs();
};

/* APIs for test cases */
class CTest {
public:
	/* create tests for specified directory */
	CTest(enum TestType, const File &, const ExeSpace & exec);
	/* deconstructor */
	~CTest();

	/* get the root directory */
	const File & get_root() const { return root; }
	/* get the space of tests */
	TestSpace & get_space() const { return *space; }

	/* load test cases from suite */
	bool load();
	/* clear the test space */
	bool clear() { space->clear(); return true; }

	/* create a new (empty) test set in the project */
	TestSet * malloc_test_set();
	/* delete the existing test set from project */
	void delete_test_set(TestSet *);
	/* whether the test set is still useful */
	bool is_accessible(TestSet *) const;

private:
	const File & root;
	enum TestType etype;
	TestSource * source;
	TestSpace * space;

	std::set<TestSet *> test_pool;
};

