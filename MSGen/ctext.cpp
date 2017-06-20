#include "ctext.h"
#include <iostream>

/** class LineReader: implement **/
LineReader::LineReader(const std::string & filename) : in(filename, std::ios::in) {
	if (!in) {
		std::cerr << "Invalid filename: \"" << filename << "\"" << std::endl;
		exit(1);
	}
}
LineReader::~LineReader() {
	in.close();
}
bool LineReader::hasNext() {
	return !(in.eof());
}
std::string LineReader::next() {
	this->roll();
	return line;
}
void LineReader::roll() {
	line = "";
	if (!in.eof()) {
		std::getline(in, line);
	}
}

/** class TextBuild: implement **/
TextBuild::TextBuild(LineReader &reader) : text(), lines() {
	lines.push_back(text.length());
	while (reader.hasNext()) {
		std::string line = reader.next();

		text += line;
		text += "\n";
		lines.push_back(text.length());
	}
}
TextBuild::~TextBuild() {}
const std::string & TextBuild::getText() const { return text; }
size_t TextBuild::numberOfCharacters() const { return text.length(); }
size_t TextBuild::numberOfLines() const { return lines.size() - 1; }
size_t TextBuild::indexOfLine(size_t line) const {
	if (line <= 0 || line >= lines.size())
		return -1;
	else return lines[line - 1];
}
size_t TextBuild::lineOfIndex(size_t index) const {
	if (index < 0 || index >= text.length()) return -1;
	else {
		// declarations
		int beg, mid, end, head, tail;
		beg = 0; end = lines.size() - 2;

		while (beg <= end) {
			mid = (beg + end) / 2;

			head = lines[mid];
			tail = lines[mid + 1];

			if (index >= head && index < tail)
				return mid + 1;
			else if (index < head)
				end = mid - 1;
			else
				beg = mid + 1;
		} /** end while binary search **/

		return -1;
	}
}

/* may be updated for LINUX */
#include <io.h>
/* derive file names in current directory */
bool derive_file_children(const std::string & file, std::vector<std::string> & list) {
	/* declarations */
	intptr_t hfile = 0L; struct _finddata_t fileinfo;
	std::string path, filename; list.clear();

	path = file + FileSeparator + "*"; bool dir = false;
	if ((hfile = _findfirst(path.c_str(), &fileinfo)) != -1) {	/* open the first child */
		dir = true;
		do {
			/* get the next child's name */
			filename = fileinfo.name;
			if (filename == ".." || filename == ".") continue;
			else filename = file + FileSeparator + filename;

			/* add to list */ list.push_back(filename);

		} while (_findnext(hfile, &fileinfo) == 0);	/* get the next child */
		_findclose(hfile); /* close the file */
	}
	return dir;
}

/* class FileItem */
void File::analysis() {
	std::vector<std::string> names;
	is_dir = derive_file_children(filename, names);

	auto beg = names.begin(), end = names.end();
	while (beg != end) {
		File & child = *(new File(*(beg++)));
		child.parent = this; children.push_back(&child);
	}
}
std::string File::to_string() {
	std::string str;

	str = "Path:\t";
	str += this->filename;
	str += "\n";

	str += "Name:\t";
	str += this->localname;
	str += "\n";

	str += "Type:\t";
	if (this->is_dir) str += "directory\n";
	else str += "file\n";

	str += "Part:\t";
	if (this->parent != nullptr)
		str += parent->filename;
	else str += "<NONE>";
	str += "\n";

	str += "CHLD:\n";
	auto beg = children.begin(), end = children.end();
	while (beg != end) {
		File * child = *(beg++);
		str += "\t";
		str += child->filename;
		str += " : ";
		str += child->localname;
		str += "\n";
	}

	return str;
}
void File::clear() {
	auto beg = children.begin(), end = children.end();
	while (beg != end) { (*beg)->clear(); delete *(beg++); }
	children.clear();
}

// APIs
bool is_space(char ch) {
	return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n';
}
bool exist_file(const std::string & path) {
	std::fstream file;
	file.open(path, std::ios::in);
	if (file) {
		file.close();
		return true;
	}
	else {
		return false;
	}
}
bool startswith(const std::string & x, const std::string & y) {
	int xn = x.length(), yn = y.length();
	int xi = 0, yi = 0;
	while (xi < xn && yi < yn) {
		if (x[xi++] != y[yi++]) return false;
	}
	return yi >= yn;
}
bool endswith(const std::string & x, const std::string & y) {
	int xn = x.length(), yn = y.length();
	while (--xn >= 0 && --yn >= 0) {
		if (x[xn] != y[yn]) return false;
	}
	return yn < 0;
}
bool trimstring(std::string & text) {
	int beg, end; char ch;

	beg = 0; end = text.length();
	while (beg < end) {
		ch = text[beg];
		if (is_space(ch)) beg++;
		else break;
	}

	while (beg < end) {
		ch = text[end - 1];
		if (is_space(ch)) end--;
		else break;
	}

	text = text.substr(beg, end - beg);
	return true;
}
int index_of(const std::string & text, char c) {
	char ch; int i = 0, n = text.length();

	while (i < n) {
		ch = text[i];
		if (ch == c) return i;
		else i++;
	}

	return -1;
}
int last_index_of(const std::string & text, char c) {
	char ch; int i = text.length() - 1;
	while (i >= 0) {
		ch = text[i];
		if (ch == c) return i;
		else i--;
	}
	return -1;
}
