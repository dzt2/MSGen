#include "ctest.h"

TestSpace::TestSpace(const CTest & m, const TestSource & src, const ExeSpace & exec)
	: project(m), source(src), test_list(), input_set(), exec_space(exec) {}
const TestCase & TestSpace::get_test(const TestCase::ID tid) const {
	if (tid >= test_list.size()) {
		std::cerr << "Invalid id: " << tid << std::endl;
		throw tid;
	}
	else {
		return *(test_list[tid]);
	}
}
bool TestSpace::create_test(const std::string & inputstr) {
	if (input_set.count(inputstr) > 0) return false;
	else {
		TestCase::ID new_id = test_list.size();
		TestCase * tc = new TestCase(*this, new_id,
			exec_space.get_exec_program(), inputstr);
		test_list.push_back(tc); input_set.insert(inputstr);
		return true;
	}

}
void TestSpace::clear() {
	auto beg = test_list.begin(), end = test_list.end();
	while (beg != end) { delete *(beg++); }
	test_list.clear(); input_set.clear();
}

TestSource::TestSource(const File & root) {
	/* initialization */
	input_dir = nullptr;
	output_dir = nullptr;
	new_output_dir = nullptr;
	suite_dir = nullptr;

	/* get children files */
	const std::vector<File *> files = root.list_files();
	auto beg = files.begin(), end = files.end();

	/* iterate all the files in the directory */
	while (beg != end) {
		/* derive the next file */
		File * file = *(beg++);
		const std::string & name
			= file->get_local_name();

		/* only when file is directory */
		if (file->is_directory()) {
			if (name == "inputs" && input_dir == nullptr) {
				input_dir = file;
			}
			else if (name == "outputs" && output_dir == nullptr) {
				output_dir = file;
			}
			else if (name == "new_outputs" && new_output_dir == nullptr) {
				new_output_dir = file;
			}
			else if (name == "suites" && suite_dir == nullptr) {
				suite_dir = file;
			}
		}

	} /* end while for root */

	  /* validation */
	std::string error = root.get_path();
	if (input_dir == nullptr) error += "/inputs/";
	else if (output_dir == nullptr) error += "/outputs/";
	else if (new_output_dir == nullptr) error += "/new_outputs/";
	else if (suite_dir == nullptr) error += "/suites/";
	else return;

	CError errors(CErrorType::InvalidArguments, "TestSource", "Undefined file: " + error);
	CErrorConsumer::consume(errors);
}

bool TestSet::has_test(TestCase::ID tid) const {
	if (tid >= space.number_of_tests()) return false;
	else return (set_vector.get_bit(tid) == BIT_1);
}
const TestCase & TestSet::get_test(TestCase::ID tid) const {
	if (tid >= space.number_of_tests() || set_vector.get_bit(tid) != BIT_1) {
		std::cerr << "Invalid tid: " << tid << std::endl;
		throw tid;
	}
	else { return space.get_test(tid); }
}
void TestSet::clear() { set_vector.clear_bytes(); number = 0; }
bool TestSet::add_test(TestCase::ID tid) {
	if (tid >= space.number_of_tests() || set_vector.get_bit(tid) == BIT_1) return false;
	else { set_vector.set_bit(tid, BIT_1); number++; return true; }
}
bool TestSet::del_test(TestCase::ID tid) {
	if (tid >= space.number_of_tests() || set_vector.get_bit(tid) == BIT_0) return false;
	else { set_vector.set_bit(tid, BIT_0); number--; return true; }
}
void TestSet::complement() {
	int length = set_vector.byte_number();
	byte * bytes = set_vector.get_bytes();
	for (int i = 0; i < length; i++)
		bytes[i] = ~(bytes[i]);
	number = space.number_of_tests() - number;
}

TestLoader::TestLoader(const TestSpace & spac) : space(spac) {
	const File & suite_dir = space.get_source().get_suite_dir();
	const std::vector<File *> & files = suite_dir.list_files();
	cur = files.begin(), end = files.end(); reader = nullptr;

	roll_next();
}
TestLoader::~TestLoader() {
	if (reader != nullptr) delete reader;
}
void TestLoader::roll_next() {
	/* initialization */
	std::string line; lineptr = nullptr;

	/* when readable and still have next file in suites */
	while (reader != nullptr || cur != end) {
		/* initialize reader to the next file */
		while (reader == nullptr && cur != end) {
			const File & file = *(*(cur++));
			reader = new LineReader(file.get_path());
		}

		/* reading lines in current suite file */
		if (reader != nullptr) {
			while (reader->hasNext()) {
				line = reader->next();
				if (line.empty()) continue;
				else { lineptr = new std::string(line); return; }
			}

			/* reading is over and release the resources */
			delete reader; reader = nullptr;
		}
	} /* end reading lines from more file */

	/* return null for no more lines */ return;
}

std::string * TestLoader_Printtokens::next_inputs() {
	/* initialization */
	std::string * lineptr = this->lineptr; roll_next();

	/* interpret the next line: (<)? path */
	if (lineptr != nullptr) {
		/* declarations */
		std::string & line = *lineptr; char ch;
		size_t k = 0, length = line.length();
		std::string path; bool is_std = false;

		/* skip the previous spaces */
		while (k < length) {
			ch = line[k];
			if (is_space(ch)) k++;
			else break;
		}
		if (k >= length) { delete lineptr; return nullptr; }

		/* '<' path */
		if (ch == '<') {
			is_std = true; k++;
			while (k < length) {
				ch = line[k];
				if (is_space(ch)) k++;
				else break;
			}
			if (k >= length) { delete lineptr; return nullptr; }
		}

		/* path */
		while (k < length) {
			ch = line[k++];
			if (is_space(ch)) break;
			else path += ch;
		}
		delete lineptr;

		/* check existence */
		if (!exist_file(space.get_source().get_input_dir().get_path() + FileSeparator + path)) return nullptr;

		path = ".." + FileSeparator + space.get_source().get_input_dir().get_local_name() + FileSeparator + path;
		if (is_std) return new std::string("< " + path + " ");
		else return new std::string(path + " ");
	}
	else return nullptr;
}
std::string * TestLoader_Replace::next_inputs() {
	/* initialization */
	std::string * lineptr = this->lineptr; roll_next();

	/* interpret as: 'xxx' 'xxx' < filepath */
	if (lineptr != nullptr) {
		/* declarations */
		std::string & line = *lineptr; char ch;
		size_t k = 0, length = line.length();
		std::string format, path, *ans = nullptr;

		/* get split */
		int e = last_index_of(line, '<');
		if (e >= 0) {
			/* get previous */
			format = line.substr(0, e);

			/* skip spaces between '<' and path */
			for (k = e + 1; k < length; k++) {
				ch = line[k];
				if (!is_space(ch)) break;
			}
			/* get path */
			path = "";
			while (k < length) {
				ch = line[k++];
				if (is_space(ch)) break;
				else path += ch;
			}

			/* validate path */
			std::string prefix = space.get_source().get_input_dir().get_path();
			prefix += FileSeparator + path;
			if (exist_file(prefix)) {
				/* construct path */
				prefix = ".." + FileSeparator;
				prefix += space.get_source().get_input_dir().get_local_name();
				path = prefix + FileSeparator + path;

				ans = new std::string(format + "< " + path + " ");
			}
		}

		/* delete a new line */
		delete lineptr; return ans;
	}
	else return nullptr;

}
std::vector<int> __int_list;
std::string * TestLoader_Schedule::next_inputs() {
	/* initialization */
	std::string * lineptr = this->lineptr; roll_next();

	/* format: int int int ... '<' path */
	if (lineptr != nullptr) {
		/* declarations */
		std::string & line = *lineptr; char ch;
		size_t k = 0, length = line.length();
		__int_list.clear(); std::string path;

		/* derive numbers before '<' */
		std::string numstr;
		while (k < length) {
			/* skip previous spaces */
			while (k < length) {
				ch = line[k];
				if (is_space(ch)) k++;
				else break;
			}
			if (k >= length) { delete lineptr; return nullptr; }
			else if (ch == '<') { k++; break; }

			/* get number string */
			numstr = "";
			while (k < length) {
				ch = line[k];
				if (is_space(ch)) break;
				else numstr += ch;
				k++;
			}

			/* get next number */
			__int_list.push_back(std::stoi(numstr));
		} /* end numbers */
		if (k >= length) { delete lineptr; return nullptr; }

		/* skip spaces between '<' and path */
		while (k < length) {
			ch = line[k];
			if (is_space(ch)) k++;
			else break;
		}
		if (k >= length) { delete lineptr; return nullptr; }

		/* derive path */
		while (k < length) {
			ch = line[k];
			if (is_space(ch)) break;
			else path += ch;
			k++;
		}
		delete lineptr;

		/* validation */
		if (!exist_file(space.get_source().get_input_dir().get_path() + FileSeparator + path)) return nullptr;
		else {
			std::string & ans = *(new std::string());
			auto beg = __int_list.begin(), end = __int_list.end();
			while (beg != end) {
				int num = *(beg++);
				ans += std::to_string(num);
				ans += " ";
			}
			ans += "< ";
			ans += ".." + FileSeparator + space.get_source().get_input_dir().get_local_name();
			ans += FileSeparator + path + " ";

			return &ans;
		}
	}
	else return nullptr;
}
std::string * TestLoader_Tcas::next_inputs() {
	/* get next line */
	std::string * lineptr = this->lineptr; roll_next();

	/* format: num num num ... num */
	if (lineptr != nullptr) {
		/* declarations */
		std::string & line = *lineptr; char ch;
		size_t k = 0, length = line.length();
		std::string numstr; __int_list.clear();

		/* derive all integers */
		while (k < length) {
			/* skip spaces */
			while (k < length) {
				ch = line[k];
				if (is_space(ch)) k++;
				else break;
			}
			if (k >= length) break;

			/* derive next numbers */
			numstr = "";
			while (k < length) {
				ch = line[k];
				if (is_space(ch)) break;
				else numstr += ch;
				k++;
			}

			/* parse integer into list */
			__int_list.push_back(std::stoi(numstr));
		}
		delete lineptr;

		/* generate input-string */
		if (!__int_list.empty()) {
			std::string & ans = *(new std::string());
			auto beg = __int_list.begin(), end = __int_list.end();
			while (beg != end) {
				int num = *(beg++);
				ans += std::to_string(num);
				if (beg != end) ans += " ";
			}
			ans += " ";
			return &ans;
		}
		/* empty line is invalid */
		else return nullptr;
	}
	else return nullptr;
}
std::string * TestLoader_Totinfo::next_inputs() {
	/* get next line */
	std::string * lineptr = this->lineptr; roll_next();

	/* format: '<' path */
	if (lineptr != nullptr) {
		/* initialization */
		std::string & line = *lineptr; char ch;
		size_t k = 0, length = line.length();
		std::string path;

		/* skip to '<' */
		while (k < length) {
			ch = line[k++];
			if (ch == '<') break;
		}
		if (k >= length) { delete lineptr; return nullptr; }

		/* skip spaces between '<' and path */
		while (k < length) {
			ch = line[k];
			if (is_space(ch)) k++;
			else break;
		}

		/* get path */
		path = "";
		while (k < length) {
			ch = line[k++];
			if (is_space(ch)) break;
			else path += ch;
		}
		delete lineptr;

		/* validation */
		if (!exist_file(space.get_source().get_input_dir().get_path() + FileSeparator + path)) return nullptr;
		else return new std::string("< .." + FileSeparator + space.get_source().get_input_dir().get_local_name() + FileSeparator + path + " ");
	}
	else return nullptr;
}

std::string * TestLoader_Flex::next_inputs() {
	/* get next line */
	std::string * lineptr = this->lineptr; roll_next();

	/* format: -P[inputstring] -F[lex.yy.c|out1] -F[error|out2] -C[V0.xxx] */
	/* target: {flex.exe} error -configlist ../inputs/filex \n cat ../exec/lex.yy.c split ../exec/error split ../exec/lex.backup */
	if (lineptr != nullptr) {
		/* get next line */
		std::string & line = *lineptr;
		int k = 0, n = line.length(); char ch;
		std::string * ans = nullptr;
		char tag; std::string content;

		/* content list */
		std::string inlist, oulist;

		/* get the inlist & oulist */
		while (k < n) {
			tag = next_tag(line, k);
			extracts(line, k, content);

			if (tag == 'P') 
				inlist = content;
			else if (tag == 'F') {
				for (int j = 0; j < content.length(); j++) {
					if (content[j] == '|') break;
					else oulist += content[j];
				}
				oulist += " split ";
			}
		}

		/* extract -o from inlist, oulist */
		extrange(inlist, oulist);

		/* generate commands */
		if (!inlist.empty()) {
			ans = new std::string();
			*ans += inlist + "\n";
			*ans += "cat " + oulist;
			*ans += " ";
		}

		delete lineptr; return ans;
	}
}
char TestLoader_Flex::next_tag(const std::string & line, int & k) {
	int n = line.length();
	while (k < n) {
		if (line[k++] == '-') break;
	}
	if (k >= n) return '\0';
	else return line[k++];
}
bool TestLoader_Flex::extracts(const std::string & line, int & k, std::string & text) {
	int n = line.length();

	text = ""; char ch;
	while (k < n) {
		if (line[k++] == '[') break;
	}
	while (k < n) {
		ch = line[k++];
		if (ch == ']') break;
		else text = text + ch;
	}

	return true;
}
bool TestLoader_Flex::extrange(std::string & inlist, std::string & oulist) {
	std::string incache;
	for (int i = 0; i < inlist.length(); i++) {
		incache += inlist[i];
		if (inlist[i] == '-') {
			incache += inlist[++i];
			if (inlist[i] == 'o') {
				while (i < inlist.length()) {
					if (inlist[i] == ' ') break;
					else i++;
				}
				incache += "./outputs "; 
				oulist += " split outputs";
			}
		}
	}
	inlist = incache;
	return true;
}

std::string * TestLoader_Gzip::next_inputs() {
	/* get next line */
	std::string * lineptr = this->lineptr; roll_next();

	/* format: -P[configs] -I[inputfile] -O[output] -X[shell]... */
	/* target: {flex.exe} error -configlist ../inputs/filex \n cat ../exec/lex.yy.c split ../exec/error split ../exec/lex.backup */
	if (lineptr != nullptr) {
		/* get next line */
		std::string & line = *lineptr;
		int k = 0, n = line.length(); char ch;
		std::string * ans = nullptr;
		char tag; std::string content;

		/* arguments */
		std::string inputstr;

		while (k < n) {
			tag = get_tag(line, k);
			get_content(line, k, content);

			if (tag == 'P') {
				inputstr += content + " ";
			}
			else if (tag == 'I') {
				inputstr += "../inputs/" + content + " ";
			}
		}

		if (content.empty()) return nullptr;
		else {
			if (is_test(inputstr)) 
				inputstr += " 2";
			else inputstr += " -c ";
			return new std::string(inputstr);
		}
	}
	else return nullptr;
}
char TestLoader_Gzip::get_tag(const std::string & line, int & k) {
	int n = line.length();

	while (k < n) {
		if (line[k++] == '-') {
			if (k < n) return line[k++];
		}
	}
	return '\0';
}
bool TestLoader_Gzip::get_content(const std::string & line, int & k, std::string & content) {
	int n = line.length();
	char ch; content = "";
	
	while (k < n) {
		if (line[k++] == '[') break;
	}
	while (k < n) {
		ch = line[k++];
		if (ch == ']') break;
		else content += ch;
	}
	
	return !content.empty();
}
bool TestLoader_Gzip::is_test(const std::string & inputstr) {
	for (int k = 0; k < inputstr.length(); k++) {
		if (inputstr[k] == '-') {
			std::string tag = "";
			while (k < inputstr.length()) {
				if (inputstr[k] == ' ') break;
				else tag += inputstr[k++];
			}

			if (tag == "-t" || 
				tag == "--test") 
				return true;
		}
	}
	return false;
}

CTest::CTest(enum TestType type, const File & dir, const ExeSpace & exec) : root(dir), etype(type), test_pool() {
	source = new TestSource(root);
	space = new TestSpace(*this, *source, exec);
}
CTest::~CTest() {
	delete space;
	delete source;

	auto beg = test_pool.begin();
	auto end = test_pool.end();
	while (beg != end)
		delete *(beg++);
	test_pool.clear();
}
bool CTest::load() {
	/* get loader */
	TestLoader *loader;
	switch (etype) {
	case print_tokens: case print_tokens2:
		loader = new TestLoader_Printtokens(*space); break;
	case replace:
		loader = new TestLoader_Replace(*space); break;
	case schedule: case schedule2:
		loader = new TestLoader_Schedule(*space); break;
	case tcas:
		loader = new TestLoader_Tcas(*space); break;
	case tot_info:
		loader = new TestLoader_Totinfo(*space); break;
	case flex:
		loader = new TestLoader_Flex(*space); break;
	case gzip:
		loader = new TestLoader_Gzip(*space); break;
	default:
		loader = new TestLoader(*space);
	}

	/* loading */
	space->clear();
	while (loader->has_next()) {
		std::string * input = loader->next_inputs();
		if (input != nullptr) {
			space->create_test(*input);
			delete input;
		}
	}

	/* return */
	delete loader;
	return true;
}
TestSet * CTest::malloc_test_set() {
	TestSet * tests = new TestSet(*space);
	test_pool.insert(tests);
	return tests;
}
void CTest::delete_test_set(TestSet * tests) {
	if (test_pool.count(tests) > 0) {
		test_pool.erase(tests);
		delete tests;
	}
}
bool CTest::is_accessible(TestSet * tests) const {
	return test_pool.count(tests) > 0;
}