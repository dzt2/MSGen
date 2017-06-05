#pragma once

/*
-File : cscript.h
-Arth : Lin Huan
-Date : May 16th 2017
-Clas :
[1] CScriptSource
[2] CScriptSpace
[3] CScriptFile
[4] CTraceScriptFile
[5] CScript
*/

#include "cfile.h"
#include "ctest.h"
#include "ctrace.h"
#include <map>

class CScript;
class CScriptSource;
class CScriptSpace;
class CScriptFile;
class CTraceScriptFile;

/* source of ../script/ */
class CScriptSource {
protected:
	/* construct source for ../script/ */
	CScriptSource(const File & d) : dir(d) {}
	/* deconstructor */
	~CScriptSource() {}
public:
	/* get ../script/ */
	const File & get_script_dir() const { return dir; }

	/* create and delete */
	friend class CScript;
private:
	/* ../script/ */
	const File & dir;
};
/* space for ../script/ of script files */
class CScriptSpace {
protected:
	/* construct space for ../script */
	CScriptSpace(const CScriptSource &, const ExeSpace &, const TestSpace &);
	/* deconstructor */
	~CScriptSpace() { clear_script(); }

	/* add a new script file */
	bool add_script(CScriptFile &);
	/* remove all script files from the space */
	void clear_script() { names.clear(); files.clear(); }
public:
	/* get the source of script */
	const CScriptSource & get_source() const { return source; }
	/* get the ../exec/ */
	const ExeSpace & get_exec() const { return exec; }
	/* get the test space of this script */
	const TestSpace & get_test() const { return test; }

	/* get the local names of script files in ../script/ */
	const std::set<std::string> & get_script_names() const { return names; }
	/* whether script file exists in the space */
	bool has_script(const std::string &) const;
	/* get the script file (for generation) in ../script/ space */
	CScriptFile & get_script(const std::string &) const;

	/* create and delete and add script file */
	friend class CScript;
private:
	/* source for ../script/ */
	const CScriptSource & source;
	/* ../exec */
	const ExeSpace & exec;
	/* test space for testing */
	const TestSpace & test;
	/* names of script files in ../script */
	std::set<std::string> names;
	/* map from script name to script file */
	std::map<std::string, CScriptFile *> files;
};
/* file for test script in ../script/ */
class CScriptFile {
protected:
	/* construct a general script file for pure testing */
	CScriptFile(const CScriptSpace & spac, const std::string & nam)
		: space(spac), name(nam), timeout(0L), echo(false) {}
	/* deconstructor */
	virtual ~CScriptFile() {}
public:
	/* get the space where the script is defined */
	const CScriptSpace & get_space() const { return space; }
	/* get the name of this script in the space */
	const std::string & get_name() const { return name; }

	/* get the time-out, 0 represent no timeout is established */
	long get_timeout() const { return timeout; }
	/* whether echo is used in script file */
	bool is_echo() const { return echo; }
	/* set time-out */
	void set_timeout(long t) { timeout = t; }
	/* set whether test information is printed to stdout */
	void set_echo(bool e) { echo = e; }

	/* generate ../script/{name}.sh based on its template */
	virtual void generate() const;

	/* create and delete */
	friend class CScript;
protected:
	const CScriptSpace & space;
	const std::string name;
	long timeout; bool echo;
};
/* final class for gettraces.sh */
class CTraceScriptFile : public CScriptFile {
protected:
	CTraceScriptFile(const CScriptSpace & spac, const std::string & nam, const CodeSpace & cs,
		const CoverageSpace & cov) : CScriptFile(spac, nam), cspace(cs), tspace(cov) {}
	~CTraceScriptFile() { clear_code_files(); }
public:
	/* get the space of code-files: ../source/ */
	const CodeSpace & get_code_space() const { return cspace; }
	/* get the space of coverage: ../traces/ */
	const CoverageSpace & get_trace_space() const { return tspace; }

	/* add a new code file in the space to the trace-script for recording its xxx.gcov in ../trace/xxx.c/ */
	bool add_code_file(const CodeFile &);
	/* clear all code files of which coverage is going to be recorded */
	void clear_code_files() { code_files.clear(); }

	/* generate ../script/{name}.sh */
	void generate() const;

	/* create and delete */
	friend class CScript;
protected:
	/* space for ../source/ */
	const CodeSpace & cspace;
	/* space for ../traces/ */
	const CoverageSpace & tspace;
	/* list of code files */
	std::set<const CodeFile *> code_files;
};
/* project for ../script */
class CScript {
public:
	/* construct script project for ../script */
	CScript(const File &, const ExeSpace &, const TestSpace &);
	/* deconstructor */
	~CScript() { clear_scripts(); delete space; delete source; }

	/* get ../ */
	const File & get_root() const { return root; }
	/* get space of ../script/ */
	const CScriptSpace & get_space() const { return *space; }

	/* create a script file in ../script */
	bool create_script(const std::string &);
	/* create a gettrace file in ../script */
	bool create_traces(const std::string &, const CodeSpace &, const CoverageSpace &);
	/* remove all script files in the space */
	bool clear_scripts();
private:
	/* ../ */
	const File & root;
	/* source of script space */
	CScriptSource * source;
	/* space for ../script */
	CScriptSpace  * space;
	/* pool for allocating memories for CScriptFile */
	std::set<CScriptFile *> file_pool;
};