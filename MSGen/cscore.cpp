#include "cscore.h"

bool ScoreVector::kill(TestCase::ID tid) {
	const TestSet & tests = function.get_tests();
	if (tests.has_test(tid)) {
		BitSeq::size_t index = function.get_index_of(tid);
		if (svec.get_bit(index) == BIT_0) {
			degree++; svec.set_bit(index, BIT_1);
			return true;
		}
		else return true;
	}
	else return false;
}

ScoreFunction::ScoreFunction(const ScoreSource & src, const TestSet & ts,
	const MutantSet & ms) : source(src), tests(ts), mutants(ms), bid_tid(), tid_bid() {
	const BitSeq & tvec = ms.get_set_vector();
	BitSeq::size_t len = tvec.bit_number();
	TestCase::ID tid = 0;

	/* construct map */
	while (tid < len) {
		if (tvec.get_bit(tid) == BIT_1) {
			tid_bid[tid] = bid_tid.size();
			bid_tid.push_back(tid);
		}
		tid = tid + 1;
	}
}
TestCase::ID ScoreFunction::get_test_id_at(BitSeq::size_t index) const {
	if (index >= bid_tid.size()) {
		CError error(CErrorType::InvalidArguments, "ScoreFunction::get_test_id_at", "Invalid index (" + std::to_string(index) + ")");
		CErrorConsumer::consume(error);
	}
	return bid_tid[index];
}
BitSeq::size_t ScoreFunction::get_index_of(TestCase::ID tid) const {
	if (tid_bid.count(tid) == 0) {
		CError error(CErrorType::InvalidArguments, "ScoreFunction::get_index_of", "Invalid test (" + std::to_string(tid) + ")");
		CErrorConsumer::consume(error);
	}
	auto iter = tid_bid.find(tid);
	return iter->second;
}

ScoreSource::ScoreSource(const CScore & p, const MutantSpace & ms, const File & rf)
	: project(p), codefile(ms.get_code_file()), mspace(ms),
	tspace(p.get_test_project().get_space()), rfile(rf), function_pool() {
	if (rfile.is_directory()) {
		CError error(CErrorType::InvalidArguments, "ScoreSource::ScoreSource", "Invalid result file: " + rfile.get_path());
		CErrorConsumer::consume(error);
	}
}
ScoreSource::~ScoreSource() {
	auto beg = function_pool.begin();
	auto end = function_pool.end();
	while (beg != end) {
		delete *(beg++);
	}
	function_pool.clear();
}
ScoreFunction * ScoreSource::create_function(
	const TestSet & tests, const MutantSet & mutants) {
	/* validation */
	const TestSpace & tspace = tests.get_space();
	const MutantSpace & mspace = mutants.get_space();
	if ((&tspace != &(this->tspace)) || (&mspace != &(this->mspace)))
		return nullptr;
	/* create a new score function */
	else {
		ScoreFunction * func = new ScoreFunction(*this, tests, mutants);
		function_pool.insert(func); return func;
	}
}
bool ScoreSource::is_accessible(ScoreFunction * func) const {
	return function_pool.count(func) > 0;
}
bool ScoreSource::delete_function(ScoreFunction * func) {
	if (function_pool.count(func) == 0) return false;
	else {
		delete func; function_pool.erase(func); return true;
	}
}

CScore::CScore(const File & root, const CMutant & mp, const CTest & tp)
	: mutproject(mp), testproject(tp), sources() {
	/* get ../score */
	dir = nullptr; auto files = root.list_files();
	auto beg = files.begin(), end = files.end();
	while (beg != end) {
		const File * file = *(beg++);
		if (file->is_directory()) {
			const std::string & name = file->get_local_name();
			if (name == "score") {
				dir = file; break;
			}
		}
	}
	if (dir == nullptr) {
		CError error(CErrorType::InvalidArguments, "CScore::CScore", "/score is not found in \"" + root.get_path() + "\"");
		CErrorConsumer::consume(error);
	}

	/* get valid code-files */
	std::map<std::string, MutantSpace *> names;
	auto mspaces = mutproject.get_spaces();
	auto mbeg = mspaces.begin(), mend = mspaces.end();
	while (mbeg != mend) {
		/* get the name (without postfix) for the next code file with mutants */
		const CodeFile * cfile = (mbeg)->first;
		MutantSpace * space = (mbeg)->second;
		std::string name = cfile->get_file().get_local_name();
		name = name.substr(0, name.length() - 2);

		/* insert the valid name of directories in ../score, which refers to mutants */
		names[name + ".txt"] = space;

		/* to the next mutant space */
		mbeg++;
	}

	/* create score function */
	files = dir->list_files();
	beg = files.begin(), end = files.end();
	while (beg != end) {
		const File * file = *(beg++);
		if (!(file->is_directory())) {
			const std::string & name = file->get_local_name();
			if (names.count(name) > 0) {
				auto iter = names.find(name);
				MutantSpace * space = iter->second;

				ScoreSource * source = new ScoreSource(*this, *space, *file);
				sources[&(space->get_code_file())] = source;
			}
		}
	}
}
CScore::~CScore() {
	auto beg = sources.begin();
	auto end = sources.end();
	while (beg != end) {
		delete (beg++)->second;
	}
	sources.clear();
}
ScoreSource & CScore::get_source(const CodeFile & cfile) const {
	if (sources.count(&cfile) == 0) {
		CError error(CErrorType::InvalidArguments, "CScore::get_source", "No mutants/results for: " + cfile.get_file().get_path());
		CErrorConsumer::consume(error);
	}
	auto iter = sources.find(&cfile);
	return *(iter->second);
}

FileScoreProducer::FileScoreProducer(const ScoreFunction & func)
	: function(func), reader(func.get_source().get_result_file().get_path()) {}
ScoreVector * FileScoreProducer::produce() {
	/* declarations */
	std::string linestr, numstr; char ch;
	TestCase::ID tid; Mutant::ID mid;
	ScoreVector * ans = nullptr; int k, n;

	/* get the next line to generate score vector */
	while (reader.hasNext()) {
		linestr = reader.next();	/* get next line */
		k = 0, n = linestr.length();

		/* get mutant id */
		numstr = "";
		while (k < n) {
			ch = linestr[k++];
			if (ch == '[') break;
			else if (ch >= '0' && ch <= '9')
				numstr += ch;
		}
		if (k >= n) continue;
		mid = std::stoul(numstr);
		ans = new ScoreVector(function, mid,
			function.get_tests().size());

		/* to the ':' */
		while (k < n) {
			ch = linestr[k++];
			if (ch == ':') break;
		}

		/* get test id from following list */
		while (k < n) {
			/* skip the 't' */
			while (k < n) {
				ch = linestr[k++];
				if (ch == 't') break;
			}
			if (k >= n) break;

			/* get test id */
			numstr = "";
			while (k < n) {
				ch = linestr[k++];
				if (is_space(ch)) break;
				else numstr += ch;
			}
			tid = std::stoul(numstr);

			/* kill by the test */
			ans->kill(tid);
		} /* end for all test id */

		break;
	}

	/* return */ return ans;
}
ScoreVector * ScoreFilter::produce() {
	ScoreVector * vec;
	while ((vec = producer.produce()) != nullptr) {
		if (vec->get_degree() == 0)			// filter equivalent 
			continue;
		else if (_template.count(vec->get_mutant()) == 0)	// filter unselected ones
			continue;
		else return vec;		// otherwise, return the vector
	}
	return nullptr;
}

CoverageVector * CoverageProducer::produce() {
	if (beg >= end) return nullptr;
	else {
		/* create coverage vector */
		CoverageVector * cvec = new CoverageVector(beg, tnum);

		/* get mutant line */
		Mutant & mutant = mspace.get_mutant(beg);
		const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & loc = mutation.get_location();
		const TextBuild * text = loc.get_file().get_text();
		if (text != nullptr) {
			size_t line = text->lineOfIndex(loc.get_bias());
			if (filecov.is_covered(line)) {	/* set coverage */
				const LineCoverage &lcov = filecov.get_line_coverage(line);
				(cvec->vec).assign(lcov.get_vector());
			}
		}

		/* iterate to next */ beg++;
		/* return */ return cvec;
	}
}
ScoreVector * CoverageScoreProducer::produce() {
	CoverageVector * cvec;
	while ((cvec = producer.produce()) != nullptr) {
		ScoreVector * svec = new ScoreVector(function, 
			cvec->get_mutant(), function.get_tests().size());
		(svec->svec).assign(cvec->get_coverage());
		consumer.consume(cvec); return svec;
	}
	return nullptr;
}


/*
int main() {
// create projects
File * dir = new File("print_tokens");
CProgram & cprog = *(new CProgram(*dir));
CTest & tprog = *(new CTest(TestType::print_tokens, *dir, cprog.get_exec()));
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
mprog.load_mutants_for(mspace, true);
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
ScoreSource & file_score = score.get_source(cfile);

// get score function
ScoreFunction & func = *(file_score.create_function(tests, mutants));
ScoreProducer producer(func); ScoreConsumer consumer(func);

// interpret
ScoreVector * vec; size_t count = 0;
while ((vec = producer.produce()) != nullptr) {
consumer.consume(vec); count++;
}
std::cout << "Producing " << count << " score vectors..." << std::endl;

// release the resources
delete &score; delete &mprog; delete &tprog; delete &cprog; delete dir;
// wait to exit
std::cout << "\nPress any key to exit......"; getchar(); return 0;
}
*/