#include "msgraph.h"
#include <fstream>

void print_MSG(const MSGraph & graph, const MuHierarchy & hierarchy, std::ostream & out) {
	const MutantSpace & space = graph.get_space();
	out << "M: " << space.number_of_mutants() << "\n";
	out << "C: " << graph.number_of_vertices() << "\n";
	out << "H: " << hierarchy.length() << "\n";

	const std::vector<MSGVertex *> & vertices = graph.get_vertices();
	auto beg = vertices.begin(), end = vertices.end(); size_t E = 0;
	while (beg != end) {
		MSGVertex & x = *(*(beg++));
		E += x.get_ou_degree();
	}

	out << "E: " << E << std::endl;
}
void stat_cluster(const MSGraph & graph, std::ostream & out) {
	const std::vector<MSGVertex *> & vertices = graph.get_vertices();
	auto beg = vertices.begin(), end = vertices.end();
	out << "Cluster\tDegree\tMutants\n";
	while (beg != end) {
		MSGVertex & vertex = *(*(beg++));
		out << vertex.get_id() << "\t"
			<< vertex.get_feature()->get_degree() << "\t"
			<< vertex.get_cluster().number_of_mutants() << "\n";
	}
	out << std::endl;
}
void stat_hierarchy(const MuHierarchy & hierarchy, std::ostream & out) {
	const std::vector<BitSeq::size_t> & deg_list = hierarchy.get_degrees();
	auto beg = deg_list.begin(), end = deg_list.end();
	out << "Degree\tClusters\tMutants\n";
	while (beg != end) {
		BitSeq::size_t deg = *(beg++); size_t mutants = 0;
		const std::vector<MSGVertex *> & level = hierarchy.get_vertices_at(deg);
		auto lbeg = level.begin(), lend = level.end();
		while (lbeg != lend) {
			MSGVertex & vertex = *(*(lbeg++));
			mutants += vertex.get_cluster().number_of_mutants();
		}
		out << deg << "\t" << level.size() << "\t" << mutants << "\n";
	}
	out << std::endl;
}

/*
int main() {
// create projects
std::string progname = "tcas";
File * dir = new File("../../SiemensSuite/" + progname);
CProgram & cprog = *(new CProgram(*dir));
CTest & tprog = *(new CTest(TestType::tcas, *dir, cprog.get_exec()));
CMutant & mprog = *(new CMutant(*dir, cprog.get_source()));
CScore & score = *(new CScore(*dir, mprog, tprog));

// load test cases from ../suites
tprog.load();
std::cout << "Load tests: " << tprog.get_space().number_of_tests() << std::endl;

// load mutants from ../muta/xxx/mschema.txt
auto mspaces = mprog.get_spaces();
auto mbeg = mspaces.begin(), mend = mspaces.end();
while (mbeg != mend) {
MutantSpace & mspace = *((mbeg++)->second);
mprog.load_mutants_for(mspace, false);
std::cout << "Load " << mspace.number_of_mutants() <<
" mutants for \"" << mspace.get_code_file().get_file().get_path() << "\"\n";
}

// get mutant set and test set
mbeg = mspaces.begin();
const CodeFile & cfile = *((mbeg)->first);
MutantSpace & mspace = mprog.get_mutants_of(cfile);
MutantSet & mutants = *(mspace.create_set());
TestSet & tests = *(tprog.malloc_test_set());
mutants.complement(); tests.complement();

// get score function
std::cout << "\n";
if (score.has_source(cfile)) {
ScoreSource & file_score = score.get_source(cfile);
ScoreFunction & func = *(file_score.create_function(tests, mutants));
ScoreProducer producer(func); ScoreConsumer consumer(func);

MSGraph graph(mspace);
MSGBuilder builder;
builder.build(graph, producer, consumer);

print_MSG(graph, builder.get_hierarchy(), std::cout);

std::ofstream cout("../../SiemensSuite/" + progname + "/clusters.txt");
stat_cluster(graph, cout); cout.close();

std::ofstream lout("../../SiemensSuite/" + progname + "/hierarchy.txt");
stat_hierarchy(builder.get_hierarchy(), lout); lout.close();
}
else {
std::cout << "No such result..." << std::endl;
}

// release the resources
delete &score; delete &mprog; delete &tprog; delete &cprog; delete dir;
// exit
std::cout << "\nPress any key to exit..." << std::endl;
getchar(); return 0;
}
*/