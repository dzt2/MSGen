#include "cscript.h"
#include <fstream>

CScriptSpace::CScriptSpace(const CScriptSource & src, const ExeSpace & exe, const TestSpace & ts)
	: source(src), exec(exe), test(ts), names(), files() {}
bool CScriptSpace::has_script(const std::string & name) const {
	return names.count(name) > 0;
}
CScriptFile & CScriptSpace::get_script(const std::string & name) const {
	if (files.count(name) == 0) {
		CError error(CErrorType::Runtime, "CScriptSpace::get_script", "Undefined name {" + name + "}");
		CErrorConsumer::consume(error);
	}
	auto iter = files.find(name);
	return *(iter->second);
}
bool CScriptSpace::add_script(CScriptFile & sfile) {
	const std::string & name = sfile.get_name();
	if (files.count(name) > 0) return false;

	files[name] = &sfile; names.insert(name);
	return true;
}

void CScriptFile::generate() const {
	/* ../exec/xxx.exe */
	std::string exec = ".." + FileSeparator;
	exec += space.get_exec().get_source().get_root().get_local_name();
	exec += FileSeparator + space.get_exec().get_source().get_exe().get_local_name();

	/* ../output/ */
	std::string outputs = ".." + FileSeparator;
	outputs += space.get_test().get_source().get_output_dir().get_local_name();
	outputs += FileSeparator;

	/* open output stream */
	std::string path = space.get_source().get_script_dir().get_path();
	path += FileSeparator + this->name; std::ofstream out(path);

	/* get the test case set */
	const TestSpace & tspace = space.get_test();
	size_t num = tspace.number_of_tests(); TestCase::ID tid = 0;

	/* iterate each test in space to generate their command */
	out << "#!/bin/sh\n";
	while (tid < num) {
		const TestCase & tc = tspace.get_test(tid++);

		if (echo)
			out << "echo \"running test " << tc.get_id() << "\"\n";

		if (timeout > 0)
			out << "timeout " << timeout << " ";

		out << exec << " " << tc.get_input_string() << " > " << outputs << tc.get_output_file() << "\n";

		// for flex script
		out << "rm lex.yy.c error lex.backup outputs\n";
	}

	/* flush and close */
	out.close();
}

bool CTraceScriptFile::add_code_file(const CodeFile & cfile) {
	const CodeSpace & cspace = cfile.get_space();
	if (&cspace != &(this->cspace)) {
		CError error(CErrorType::Runtime, "CTraceScriptFile::add_code_file", "Unmatched code file {" + cfile.get_file().get_path() + "}");
		CErrorConsumer::consume(error);
	}
	else if (code_files.count(&cfile)) return false;
	else {
		code_files.insert(&cfile); return true;
	}
}
void CTraceScriptFile::generate() const {
	/* ../exec/xxx.trace.exe */
	std::string exec = ".." + FileSeparator;
	exec += space.get_exec().get_source().get_root().get_local_name() + FileSeparator;
	exec += space.get_exec().get_source().get_trace_exe().get_local_name();

	/* get ../traces/ */
	std::string trace = ".." + FileSeparator;
	trace += tspace.get_source().get_trace_dir().get_local_name();
	trace += FileSeparator;

	/* get ../source/ */
	std::string source = ".." + FileSeparator;
	source += cspace.get_source().get_directory().get_local_name();
	source += FileSeparator;

	/* ../output/ */
	std::string outputs = ".." + FileSeparator;
	outputs += space.get_test().get_source().get_output_dir().get_local_name();
	outputs += FileSeparator;

	/* open output stream */
	std::string path = space.get_source().get_script_dir().get_path();
	path += FileSeparator + this->name; std::ofstream out(path);

	/* initialization: mkdir in ../traces/ */
	out << "#!/bin/sh\n";
	out << "cd " << trace << "\n";
	out << "rm -rf *\n";
	auto beg = code_files.begin();
	auto end = code_files.end();
	while (beg != end) {
		const CodeFile & cfile = *(*(beg++));
		const std::string & name = cfile.get_file().get_local_name();
		//out << "if [ ! -d \"" << name << "\" ]; then\n";
		//out << "\tmkdir " << name << "\nfi\n";
		out << "mkdir " << name << "\n";
	}

	/* cd ../source/ */ out << "\ncd " << source << "\n";

	/* generate echo-rm-exec-gcov(s) list */
	const TestSpace & tests = space.get_test();
	TestCase::ID tid = 0;
	size_t tnum = tests.number_of_tests();
	while (tid < tnum) {
		const TestCase & tc = tests.get_test(tid++);

		/* echo-rm-exec-gcov for each test */
		if (echo) out << "echo \"running test " << tc.get_id() << "\"\n";
		out << "rm *.gcda *.gcov\n";
		if (timeout > 0) out << "timeout " << timeout << " ";
		out << exec << " " << tc.get_input_string() << " > " << outputs << tc.get_output_file() << "\n";
		// out << "gcov *.c\n";
		// for flex
		out << "gcov flex.c\n";

		/* copy each code file's gcov to the ../trace/xxx.c/txx */
		beg = code_files.begin(), end = code_files.end();
		while (beg != end) {
			const CodeFile & cfile = *(*(beg++));
			const std::string & name = cfile.get_file().get_local_name();

			out << "if [ -f \"" << name << ".gcov\" ]; then\n";
			out << "\tcp " << name << ".gcov " << trace << name
				<< FileSeparator << tc.get_output_file() << "\n";
			out << "fi\n";
		} /* end while: iterate gcov file for each code file in the list to ../traces/ */

		// for flex script
		out << "rm lex.yy.c error lex.backup outputs\n";

		out << "\n";
	} /* end while: output all commands for every test in space */

	  /* flush and close */
	out.close();
}

CScript::CScript(const File & dir, const ExeSpace & exec, const TestSpace & tests) : root(dir) {
	if (!root.is_directory()) {
		CError error(CErrorType::InvalidArguments, "CScript::CScript", "Not directory: " + root.get_path());
		CErrorConsumer::consume(error);
	}

	/* get source for ../script */
	auto files = root.list_files(); source = nullptr;
	auto beg = files.begin(), end = files.end();
	while (beg != end) {
		const File & file = *(*(beg++));
		if (file.is_directory()) {
			if (file.get_local_name() == "script") {
				source = new CScriptSource(file);
				break;
			}
		}
	}

	/* validation source */
	if (source == nullptr) {
		CError error(CErrorType::InvalidArguments, "CScript::CScript", "script/ is not found: " + root.get_path());
		CErrorConsumer::consume(error);
	}

	space = new CScriptSpace(*source, exec, tests);
}
bool CScript::create_script(const std::string & name) {
	if (space->has_script(name)) return false;
	else {
		CScriptFile * sfile = new CScriptFile(*space, name);
		space->add_script(*sfile); file_pool.insert(sfile);
		return true;
	}
}
bool CScript::create_traces(const std::string & name,
	const CodeSpace & cspace, const CoverageSpace & tspace) {
	if (space->has_script(name)) return false;
	else {
		CTraceScriptFile * sfile = new CTraceScriptFile(*space, name, cspace, tspace);
		space->add_script(*sfile); file_pool.insert(sfile); return true;
	}
}
bool CScript::clear_scripts() {
	space->clear_script();
	auto beg = file_pool.begin(), end = file_pool.end();
	while (beg != end) {
		CScriptFile * sfile = *(beg++);
		delete sfile;
	}
	file_pool.clear();
	return true;
}

/*
int main() {
// initialization
File * dir = new File("schedule");
CProgram & program = *(new CProgram(*dir));
CTest & test = *(new CTest(TestType::schedule, *dir, program.get_exec()));
CTrace & trace = *(new CTrace(*dir, program.get_source(), test.get_space()));
CScript & script = *(new CScript(*dir, program.get_exec(), test.get_space()));

// load test suite into test cases
test.load();
script.create_script("runall.sh");
script.create_traces("gettrace.sh", program.get_source(), trace.get_space());

// generate runall.sh
const CScriptSpace & script_space = script.get_space();
CScriptFile & runall = script_space.get_script("runall.sh");
runall.set_echo(true); runall.set_timeout(3); runall.generate();

// generate gettrace.sh
CTraceScriptFile * traceptr = (CTraceScriptFile *)(&(script_space.get_script("gettrace.sh")));
CTraceScriptFile & gettrace = *traceptr;
auto cfiles = program.get_source().get_code_set();
auto beg = cfiles.begin(), end = cfiles.end();
while (beg != end) {
const CodeFile & cfile = *(*(beg++));
gettrace.add_code_file(cfile);
}
gettrace.set_echo(true);  gettrace.generate();

// release all resources
delete &script; delete &test; delete &trace; delete &program; delete dir;

// exit
std::cout << "No exception..." << std::endl;
getchar(); return 0;
}
*/