#include "cfunc.h"
#include <fstream>
#include <algorithm>

const CodeLocation & CFunction::get_definition_point() const {
	if (define_point == nullptr) {
		CError error(CErrorType::InvalidArguments,
			"CFunction::get_definition_point()",
			func_name + " is not defined");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else return *define_point;
}

const std::set<const CFunction *> & CFunctionGraph::get_callers_of(const CFunction & x) const {
	if (callers.count(&x) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"CFunctionGraph::get_callers_of(x)", 
			"no callers for " + x.get_name());
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = callers.find(&x);
		return *(iter->second);
	}
}
const std::set<const CFunction *> & CFunctionGraph::get_callees_of(const CFunction & x) const {
	if (callees.count(&x) == 0) {
		CError error(CErrorType::InvalidArguments,
			"CFunctionGraph::get_callee_of(x)",
			"no callees for " + x.get_name());
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = callees.find(&x);
		return *(iter->second);
	}
}
void CFunctionGraph::add_calling(const CFunction & x, const CFunction & y) {
	if (callees.count(&x) == 0)
		callees[&x] = new std::set<const CFunction *>();
	if (callers.count(&y) == 0)
		callers[&y] = new std::set<const CFunction *>();

	auto xiter = callees.find(&x);
	(xiter->second)->insert(&y);
	auto yiter = callers.find(&y);
	(yiter->second)->insert(&x);
}
void CFunctionGraph::clear_graph() {
	auto beg1 = callees.begin();
	auto end1 = callees.end();
	while (beg1 != end1) 
		delete (beg1++)->second;
	callees.clear();

	auto beg2 = callers.begin();
	auto end2 = callers.end();
	while (beg2 != end2)
		delete (beg2++)->second;
	callers.clear();

}

const CFunction & CFunctionSpace::get_function(const std::string & name) const {
	if (functions.count(name) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"CFunctionSpace::get_function(name)", 
			name + " is not defined or declared");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = functions.find(name);
		return *(iter->second);
	}
}
void CFunctionSpace::add_function(const std::string & name) {
	if (functions.count(name) > 0) {
		CError error(CErrorType::InvalidArguments,
			"CFunctionSpace::get_function(name)",
			name + " has been defined or declared");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		CFunction * func = new CFunction(*this, name);
		functions[name] = func; func_names.insert(name);
	}
}
void CFunctionSpace::add_function(const std::string & name, size_t begin, size_t length) {
	if (functions.count(name) > 0) {
		CError error(CErrorType::InvalidArguments,
			"CFunctionSpace::get_function(name)",
			name + " has been defined or declared");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		CodeLocation loc(file, begin, length);
		CFunction * func = new CFunction(*this, name, loc);
		functions[name] = func; func_names.insert(name);
	}
}
void CFunctionSpace::clear_functions_and_graph() {
	/* delete functions */
	auto beg = functions.begin();
	auto end = functions.end();
	while (beg != end)
		delete (beg++)->second;

	/* clear sets */
	functions.clear();
	func_names.clear();
	graph->clear_graph();
}
bool CFunctionSpace::find_function_at(size_t bias, std::string & name) const {
	auto beg = functions.begin(), end = functions.end();
	while (beg != end) {
		CFunction * func = beg->second;
		if (func->is_defined()) {
			const CodeLocation & loc = func->get_definition_point();
			if (bias >= loc.get_bias() && 
				bias <= loc.get_bias() + loc.get_length()) {
				name = beg->first; return true;
			}
		}
		beg++;
	}
	return false;
}

void CFunctionLoader::load_definition(const File & fun, CFunctionSpace & space) {
	// declarations
	std::string line, name, beg;
	const CodeFile & cfile = space.get_code_file();
	std::vector<size_t> begins;
	std::map<size_t, std::string> names;

	/* collect information */
	std::ifstream in(fun.get_path());
	while (std::getline(in, line)) {
		int k = 0, n = line.length(); 

		// skip spaces 
		while (k < n) {
			if (isspace(line[k])) k++;
			else break;
		}
		if (k >= n) continue;

		// derive name
		name = "";
		while (k < n) {
			if (isspace(line[k])) break;
			else name += line[k++];
		}
		
		// skip begin 
		beg = "";
		while (k < n) {
			if (isspace(line[k])) k++;
			else break;
		}
		while (k < n) {
			if (isspace(line[k])) break;
			else beg += line[k++];
		}

		// skip name length 
		beg = "";
		while (k < n) {
			if (isspace(line[k])) k++;
			else break;
		}
		while (k < n) {
			if (isspace(line[k])) break;
			else beg += line[k++];
		}
		
		// derive line_number
		beg = "";
		while (k < n) {
			if (isspace(line[k])) k++;
			else break;
		}
		while (k < n) {
			if (isspace(line[k])) break;
			else beg += line[k++];
		}

		// invalid function line 
		if (name.empty()) {
			CError error(CErrorType::InvalidArguments,
				"CFunctionLoader::load_definition()",
				"Unable to explain: \"" + line + "\"");
			CErrorConsumer::consume(error); 
			exit(CErrorType::InvalidArguments);
		}
		size_t line_num = std::atoi(beg.c_str());
		size_t line_beg = cfile.get_text()->indexOfLine(line_num);

		/* record line and function */
		begins.push_back(line_beg); names[line_beg] = name;
	}
	in.close();

	// resort the line-begins 
	begins.push_back(cfile.get_text()->numberOfCharacters());
	std::sort(begins.begin(), begins.end());
	// create function definitions 
	space.clear_functions_and_graph();
	for (int i = 0; i < begins.size() - 1; i++) {
		size_t beg = begins[i];
		size_t end = begins[i + 1];
		auto iter = names.find(beg);
		const std::string & name = iter->second;
		space.add_function(name, beg, end - beg);
	}

	
}
void CFunctionLoader::load_call_graph(const File & cgr, CFunctionGraph & graph) {
	// declarations 
	std::string line, text, token;
	std::ifstream in(cgr.get_path());
	CFunctionSpace & space = graph.get_space();

	// get text 
	graph.clear_graph();
	while (std::getline(in, line))
		text += line + "\n";

	std::string head = "";
	int k = 0, n = text.length();
	while (k < n) {
		// skip spaces
		while (k < n) {
			if (isspace(text[k])) k++;
			else break;
		}
		
		// non-end-of-file
		if (k < n) {
			// get token 
			token = "";
			while (k < n) {
				if (isspace(text[k])) break;
				else token += text[k++];
			}

			// update head for caller
			if (token[0] == '@') head = token.substr(1);
			// insert new function call relation
			else {
				const CFunction & caller = space.get_function(head);
				if (!space.has_function(token))
					space.add_function(token);
				const CFunction & callee = space.get_function(token);

				graph.add_calling(caller, callee);
			}
		}
	} // end while
}

CFuncProject::CFuncProject(const File & dir) 
	: root(nullptr), spaces(), loader() {
	const std::vector<File *> & files = dir.list_files();
	for (int i = 0; i < files.size(); i++) {
		File & child = *files[i];
		if (child.get_local_name() == "muta") {
			root = &child; break;
		}
	}
	if (root == nullptr) {
		CError error(CErrorType::InvalidArguments, 
			"CFuncProject::CFuncProject", "../muta is not found");
		CErrorConsumer::consume(error); 
		exit(CErrorType::InvalidArguments);
	}
}
CFunctionSpace & CFuncProject::get_function_space(const CodeFile & file) const {
	if (spaces.count(&file) == 0) {
		CError error(CErrorType::InvalidArguments,
			"CFuncProject::get_function_space", 
			file.get_file().get_local_name() 
			+ " does not load for its space");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = spaces.find(&file); return *(iter->second);
	}
}
void CFuncProject::clear_spaces() {
	auto beg = spaces.begin(), end = spaces.end();
	while (beg != end) delete (beg++)->second;
	spaces.clear();
}
bool CFuncProject::load_functions_for(const CodeFile & file) {
	// invlaid code-file
	if (file.get_type() != CodeFileType::Source
		&& file.get_type() != CodeFileType::Preprocessed) {
		CError error(CErrorType::InvalidArguments,
			"CFuncProject::load_functions_for()",
			file.get_file().get_local_name()
			+ " is not valid code-file");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}

	// get code-file name
	std::string filename = file.get_file().get_local_name();
	filename = filename.substr(0, filename.length() - 2);

	// find ../muta/xxx target directory
	File * tdir = nullptr;
	const std::vector<File *> & files = root->list_files();
	for (int i = 0; i < files.size(); i++) {
		File & child = *(files[i]);
		if (child.get_local_name() == filename) {
			tdir = &child; break;
		}
	}
	if (tdir == nullptr) {
		CError error(CErrorType::InvalidArguments,
			"CFuncProject::load_functions_for()",
			filename + " is not found in ../muta/");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}

	// get ../muta/__xxx.fun and __xxx.cgr
	File * fun = nullptr; File * cgr = nullptr;
	const std::vector<File *> & tfiles = tdir->list_files();
	for (int i = 0; i < tfiles.size(); i++) {
		File * child = tfiles[i];
		if (child->get_local_name() == "__" + filename + ".fun") {
			fun = child;
		}
		else if (child->get_local_name() == "__" + filename + ".cgr") {
			cgr = child;
		}
	}

	// load functions
	if (fun != nullptr) {
		// create a new space
		if (spaces.count(&file) == 0) {
			CFunctionSpace * space = new CFunctionSpace(*this, file);
			spaces[&file] = space;
		}

		// load information 
		auto iter = spaces.find(&file);
		CFunctionSpace & space = *(iter->second);
		loader.load_definition(*fun, space);
		if (cgr != nullptr)
			loader.load_call_graph(*cgr, space.get_function_call_graph());

		return true;	// load end
	}
	else return false;
}
