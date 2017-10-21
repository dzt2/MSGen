#pragma once

/*
	-file: cfunc.h
	-purp: to define interfaces to access functions in code-file
	-arth: Lin Huan
	-date: oct 21th 2017
	-clas:
		CFunction
		CFunctionGraph
		CFunctionSpace
		CFunctionLoader
		CFunctionManage
*/

#include "cfile.h"
#include <map>

class CFunction;
class CFunctionGraph;
class CFunctionSpace;
class CFunctionLoader;
class CFuncProject;

/* function is module in code-file (defined or declared for used) */
class CFunction {
protected:
	/* create a function (not defined in space) */
	CFunction(CFunctionSpace & spac, const std::string & name) 
		: space(spac), func_name(name), define_point(nullptr) {}
	CFunction(CFunctionSpace & spac, const std::string & name, const CodeLocation & loc)
		: space(spac), func_name(name) {
		define_point = new CodeLocation(loc);
	}
	/* deconstructor */
	~CFunction() { if (define_point != nullptr) delete define_point; }

public:
	/* get the space where the function is defined */
	CFunctionSpace & get_space() const { return space; }
	/* get the function name */
	const std::string & get_name() const { return func_name; }
	/* whether this function is defined */
	bool is_defined() const { return define_point != nullptr; }
	/* get the location where the function is defined */
	const CodeLocation & get_definition_point() const;

	/* new & delete */
	friend class CFunctionSpace;

private:
	CFunctionSpace & space;
	const std::string func_name;
	CodeLocation * define_point;
};
/* function call graph */
class CFunctionGraph {
protected:
	/* create a function call graph */
	CFunctionGraph(CFunctionSpace & spac)
		: space(spac), callees(), callers() {}
	/* deconstructor */
	~CFunctionGraph() { clear_graph(); }

	/* clear the graph */
	void clear_graph();
	/* add a new call-relation from x to y: x calls y */
	void add_calling(const CFunction &, const CFunction &);

public:
	/* get the space of the function call graph */
	CFunctionSpace & get_space() const { return space; }
	/* whether there are functions that call this one */
	bool has_callers_of(const CFunction & x) const { return callers.count(&x) > 0; }
	/* whether there are functions that are called by the one */
	bool has_callees_of(const CFunction & x) const { return callees.count(&x) > 0; }
	/* get the functions that call the target one */
	const std::set<const CFunction *> & get_callers_of(const CFunction &) const;
	/* get the functions that are called by target one  */
	const std::set<const CFunction *> & get_callees_of(const CFunction &) const;

	/* new & delete */
	friend class CFunctionSpace;
	/* add-calling or clear-graph */
	friend class CFunctionLoader;

private:
	/* function space to provide function */
	CFunctionSpace & space;
	/* function to those that are called by it */
	std::map<const CFunction *, std::set<const CFunction *> *> callees;
	/* function to those that call it */
	std::map<const CFunction *, std::set<const CFunction *> *> callers;
};
/* function space for a specific code file */
class CFunctionSpace {
protected:
	/* create a function space for specified code file */
	CFunctionSpace(const CFuncProject & pro, const CodeFile & fil)
		: project(pro), file(fil), func_names(), functions() {
		graph = new CFunctionGraph(*this);
	}
	/* deconstructor */
	~CFunctionSpace() { clear_functions_and_graph(); delete graph; }

	/* add a undefined function in space */
	void add_function(const std::string &);
	/* add a defined function in the space */
	void add_function(const std::string &, size_t begin, size_t length);
	/* clear all the functions and call-graph */
	void clear_functions_and_graph();

public:
	/* get the project where this space is defined */
	const CFuncProject & get_project() const { return project; }
	/* get the source file where these functions are declared or defined */
	const CodeFile & get_code_file() const { return file; }

	/* get the set of functions names in the space */
	const std::set<std::string> & get_function_names() const { return func_names; }
	/* whether there is function to the name declared or defined in the space */
	bool has_function(const std::string & name) const { return func_names.count(name); }
	/* get the function from space by its name */
	const CFunction & get_function(const std::string &) const;
	/* get the function-call-graph */
	CFunctionGraph & get_function_call_graph() const { return *graph; }

	/* get the number of functions */
	size_t size() const { return functions.size(); }
	/* find the name of function where the byte is located, return false if finding fails */
	bool find_function_at(size_t bias, std::string & name) const;

	/* new & delete */
	friend class CFuncProject;
	/* add-function & clear_functions_and_graph */
	friend class CFunctionLoader;

private:
	/* project where the space is defined */
	const CFuncProject & project;
	/* code file where functions are declared */
	const CodeFile & file;
	/* set of function names */
	std::set<std::string> func_names;
	/* map from name to function object */
	std::map<std::string, CFunction *> functions;
	/* function-call-graph */
	CFunctionGraph * graph;
};
/* loader for functions */
class CFunctionLoader {
protected:
	/* constructor */
	CFunctionLoader() {}
	/* deconstructor */
	~CFunctionLoader() {}

public:
	/* load definitions from ../muta/xxx/xxx.fun */
	void load_definition(const File &, CFunctionSpace &);
	/* load function-call-graph from ../muta/xxx/xxx.cgr */
	void load_call_graph(const File &, CFunctionGraph &);

	/* new & delete */
	friend class CFuncProject;

private:
};

/* function project */
class CFuncProject {
public:
	/* create a function project */
	CFuncProject(const File &);
	/* deconstrucotr  */
	~CFuncProject() { clear_spaces(); }

	/* whether the code-file has been loaded */
	bool has_functions_of(const CodeFile & file) const { return spaces.count(&file) > 0; }
	/* load the function space for the code file */
	bool load_functions_for(const CodeFile &);
	/* get the function space of the code file (error when not loaded) */
	CFunctionSpace & get_function_space(const CodeFile &) const;

protected:
	/* clear all the function space from project */
	void clear_spaces();

private:
	/* director for ../muta/ */
	const File * root;
	/* map from file to corresponding funciton space */
	std::map<const CodeFile *, CFunctionSpace *> spaces;
	/* to load function-space */
	CFunctionLoader loader;
};
