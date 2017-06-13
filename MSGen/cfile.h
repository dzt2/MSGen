#pragma once

/*
-File : cfile.h
-Arth : Lin Huan
-Date : May 15th 2017
-Purp : provide interfaces to access files of source code and executionable programs.
-Clas :
[1] class ExeSource;
[2] class ExeSpace;
[3] class ExeFile;

[4] class CodeSource;
[5] class CodeSpace;
[6] class CodeFile;
enum CodeFileType;

[7] class CProgram
[8] class CodeLoader
*/

#include "ctext.h"
#include "cerror.h"
#include <set>

// declarations
class ExeSource;
class ExeSpace;
class ExeFile;

class CodeSource;
class CodeSpace;
class CodeFile;
enum CodeFileType;

class CodeLoader;
class CProgram;

/* source for ../exec/ */
class ExeSource {
protected:
	/* create exec-source by its directory */
	ExeSource(const File &);
	/* deconstructor */
	~ExeSource() {}
public:
	/* get directory for ../exec/ */
	const File & get_root() const { return root; }
	/* get file of ../exec/xxx.exe */
	const File & get_exe() const { return *program; }
	/* get file of ../exec/xxx.trace.exe */
	const File & get_trace_exe() const { return *trace; }

	/* create and delete */
	friend class CProgram;
private:
	/* ../exec */
	const File & root;
	/* ../exec/xxx.exe */
	const File * program;
	/* ../exec/xxx.trace.exe */
	const File * trace;
};
/* space for executionable programs */
class ExeSpace {
protected:
	/* create a new exec-space based on its source */
	ExeSpace(const CProgram &, const ExeSource & src);
	/* deconstructor */
	~ExeSpace() {}

public:

	const CProgram & get_project() const { return project; }
	/* get the source of ../exec/ */
	const ExeSource & get_source() const { return source; }
	/* get program of ../exec/xxx.c */
	const ExeFile & get_exec_program() const { return *exec_program; }
	/* get program of ../exec/xxx.trace.exe */
	const ExeFile & get_trace_program() const { return *trace_program; }

	friend class CProgram;
private:
	/* program of this exec space */
	const CProgram & project;
	/* ../exec/ */
	const ExeSource & source;
	/* ../exec/xxx.exe */
	ExeFile * exec_program;
	/* ../exec/xxx.trace.exe */
	ExeFile * trace_program;
};
/* file for executionable programs */
class ExeFile {
protected:
	/* create a exe-file from xxx.exe */
	ExeFile(const ExeSpace &, const File & exe);
	/* deconstructor */
	~ExeFile() {}
public:
	/* get the space of this exe-file */
	const ExeSpace & get_space() const { return space; }
	/* get ../exec/xxx.exe */
	const File & get_file() const { return file; }
	/* create and delete */
	friend class ExeSpace;
private:
	/* space for ../exec/ */
	const ExeSpace & space;
	/* ../exec/xxx.exe */
	const File & file;
};

/* source for ../source(.orig)/ */
class CodeSource {
protected:
	/* construct code source in directory */
	CodeSource(const File &);
	/* deconstructor */
	~CodeSource() {}

public:
	/* get ../source(.orig)/ */
	const File & get_directory() const { return dir; }

	friend class CProgram;
private:
	/* ../source(.orig)/ */
	const File & dir;
};
/* code file */
class CodeFile {
protected:
	/* construct a new code-file with  */
	CodeFile(const CodeSpace &, const File &);
	/* deconstructor */
	~CodeFile() { if (text != nullptr) delete text; }

	/* text to the code of this file */
	const TextBuild * text;
public:
	/* get the type of this code file */
	enum CodeFileType get_type() const { return type; }
	/* get the file of this code */
	const File & get_file() const { return source; }
	/* get the space where the code file is defined */
	const CodeSpace & get_space() const { return space; }
	/* get the text of the source code */
	const TextBuild * get_text() const { return text; }

	/* get code space */
	friend class CodeSpace;
	/* set text */
	friend class CodeLoader;
private:
	/* code space where this file is defined */
	const CodeSpace & space;
	/* file of this code */
	const File & source;
	/* type of this code file */
	enum CodeFileType type;
};
/* type for code file */
enum CodeFileType {
	Header,
	Source,
	Preprocessed,
};
/* loader to load text in code file */
class CodeLoader {
protected:
	/* code loader */
	CodeLoader() {}
	/* deconstructor */
	~CodeLoader() {}

	/* load text for the code file */
	void load_in(CodeFile & cfile) const;
	/* release the text from memory for code file */
	void release(CodeFile & cfile) const;

	/* create and delete */
	friend class CodeSpace;
};
/* space for code files */
class CodeSpace {
protected:
	/* construct a code space based on its source */
	CodeSpace(const CProgram &, const CodeSource &);
	/* deconstructor */
	~CodeSpace();
public:
	/* get the project of this code space */
	const CProgram & get_project() const { return project; }
	/* get the source of the space */
	const CodeSource & get_source() const { return source; }
	/* get the set of source code files in the space */
	const std::set<CodeFile *> & get_code_set() const { return code_set; }

	/* whether specified code file has been loaded */
	bool is_loaded(const CodeFile & cfile) const {
		return cfile.get_text() != nullptr;
	}
	/* load text into code file */
	void load(CodeFile & cfile) const {
		if (cfile.get_text() == nullptr)
			loader.load_in(cfile);
	}
	/* remove the text from memory for code file */
	void release(CodeFile & cfile) const {
		if (cfile.get_text() != nullptr)
			loader.release(cfile);
	}

	friend class CProgram;
private:
	/* project for the code space */
	const CProgram & project;
	/* source of code space */
	const CodeSource & source;
	/* set of code files */
	std::set<CodeFile *> code_set;
	/* for loading text in code file */
	const CodeLoader loader;
};

/* program for both exe and code */
class CProgram {
public:
	/* construct program space in the given siemens suite project */
	CProgram(const File &);
	/* deconstructor */
	~CProgram();

	/* get ../ */
	const File & get_root() const { return root; }
	/* get ../exec/ */
	const ExeSpace & get_exec() const { return *exe_space; }
	/* get ../source/ */
	const CodeSpace & get_source() const { return *source; }
	/* get ../source.orig/ */
	const CodeSpace & get_source_orig() const { return *source_orig; }

private:
	const File & root;
	ExeSource * exe_source;
	ExeSpace * exe_space;
	CodeSource * code_source;
	CodeSpace * source;
	CodeSource * orig_source;
	CodeSpace * source_orig;
};

/* location referring to substring in code file */
class CodeLocation {
public:
	/* construct a location  */
	CodeLocation(const CodeFile & cfile, size_t b, size_t n)
		: file(cfile), bias(b), length(n) {}
	/* deconstructor */
	~CodeLocation() {}

	/* get the index to first character of the location */
	size_t get_bias() const { return bias; }
	/* get the number of characters in the location */
	size_t get_length() const { return length; }
	/* get the file referred by this code file */
	const CodeFile & get_file() const { return file; }
	/* get the code text at specified */
	std::string get_text_at() const;

private:
	/* file referred by this location */
	const CodeFile & file;
	size_t bias, length;
};

