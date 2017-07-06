#include "ctest.h"
#include "ctrace.h"
#include "cscript.h"

/*
int main() {
// create projects
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string name = "tot_info"; 
	File & dir = *(new File(prefix + name)); CProgram & program = *(new CProgram(dir));
	CTest & tproject = *(new CTest(TestType::tot_info, dir, program.get_exec()));
	CTrace & cproject = *(new CTrace(dir, program.get_source(), tproject.get_space()));
	
	// load tests from ../suites/
	tproject.load();
	const TestSpace & tspace = tproject.get_space();
	std::cout << "Generate " << tspace.number_of_tests() << " tests.\n";
	
	// generate script file
	CScript & script = *(new CScript(dir, program.get_exec(), tspace));
	script.create_script("runall.sh");
	script.create_traces("gettrace.sh", program.get_source(), cproject.get_space());
	CScriptFile & runall = script.get_space().get_script("runall.sh");
	CTraceScriptFile & gettrace = *((CTraceScriptFile *)(&(script.get_space().get_script("gettrace.sh"))));

	// establish arguments
	runall.set_echo(true); runall.set_timeout(1);
	gettrace.set_echo(true); gettrace.set_timeout(0);
	auto cfiles = program.get_source().get_code_set();
	auto beg = cfiles.begin(), end = cfiles.end();
	while (beg != end) {
		gettrace.add_code_file(*(*(beg++)));
	}

	// generate script content
	runall.generate(); gettrace.generate();

	// release projects
	delete &script; delete &cproject;
	delete &tproject; delete &program; delete &dir;
	std::cout << "\nPress any key to exit..." << std::endl;
	getchar(); return 0;
}
*/