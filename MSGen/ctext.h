#pragma once
/**
* File: text.h
* 	-Aim: to provide interfaces for access to file text
* 	-Cls:
* 		[1] class LineReader;
* 		[2] class TextBuild;
* 	-Dat: March 10th, 2017
* 	-Art: Lin Huan
* **/

#include <string>
#include <vector>
#include <fstream>

class LineReader;
class TextBuild;
class File;

/* string analysis methods */
/* is space character */
extern bool is_space(char);
/* whether file exist */
extern bool exist_file(const std::string &);
/* whether first string begins with another */
extern bool startswith(const std::string &, const std::string &);
/* whether first string ends with another */
extern bool endswith(const std::string &, const std::string &);
/* remove spaces in both head and tail of the string */
extern bool trimstring(std::string &);
/* get the index of first occurrence of character in string, if no occurrence, return -1 */
extern int index_of(const std::string &, char);
/* get the index of last occurrence of character in string, if no occurrence, return -1 */
extern int last_index_of(const std::string &, char);

/* derive file names in current directory and return whether file is directory */
extern bool derive_file_children(const std::string & file, std::vector<std::string> & list);
/* separator between file path */
static const std::string FileSeparator = "/";

/*
* Reader to retrieve text from file line by line
* */
class LineReader {
public:
	/* constructor */
	LineReader(const std::string &);
	/* deconstructor */
	~LineReader();

	/* whether there has next line */
	bool hasNext();
	/* get next line (dynamically allocate) */
	std::string next();

private:
	/* input stream for file */
	std::ifstream in;
	/* pointer to the next line */
	std::string line;

	/* update pointer to the next line */
	void roll();
};
/*
* Build for text to retrieve character by their lines
* */
class TextBuild {
public:
	/* constructor */
	TextBuild(LineReader &);
	/* deconstructor */
	~TextBuild();

	/* number of characters in text */
	size_t numberOfCharacters() const;
	/* number of lines in text */
	size_t numberOfLines() const;
	/* index to the first character in specific line (start from 1) */
	size_t indexOfLine(size_t line) const;
	/* line in which the character of index is located (index start from 0) */
	size_t lineOfIndex(size_t index) const;
	/* get the text retrieved from file */
	const std::string & getText() const;
private:
	/* original text from source code */
	std::string text;
	/* list of indexes to the */
	std::vector<std::size_t> lines;

};
/* object for File */
class File {
public:
	/* construct file for specified path */
	File(const std::string & name) : filename(name), parent(nullptr), children() {
		analysis();

		int index = filename.length(); localname = filename;
		while (--index >= 0)
			if (filename[index] == '/' || filename[index] == '\\')
				break;
		index++;

		localname = filename.substr(index, filename.length() - index);
	}
	/* remove local file items */
	~File() { clear(); }

	/* get the parent of this file */
	const File * get_parent() const { return parent; }
	/* get the path to this file  */
	const std::string & get_path() const { return filename; }
	/* get the local file name */
	const std::string & get_local_name() const { return localname; }
	/* is this file a directory */
	bool is_directory() const { return is_dir; }
	/* whether this file exists */
	bool exist() const { return is_dir || exist_file(filename); }
	/* files in current directory */
	const std::vector<File *> & list_files() const { return children; }

	/* update the list of files under the directory */
	void update_list() { clear(); analysis(); }

	/* print information of this file */
	std::string to_string();
private:
	/* file path */
	std::string filename;
	/* local program name */
	std::string localname;
	/* parent file */
	const File * parent;
	/* is this file a directory */
	bool is_dir;
	/* sub-files in current directory */
	std::vector<File *> children;

	/* check out whether this file is directory and */
	void analysis();
	/* clear all file-images in the directory */
	void clear();
};
