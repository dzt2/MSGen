#include "ctrace.h"

void LineCoverage::cover(const BitSeq::size_t tid) {
	vec.set_bit(tid, BIT_1);
}

FileCoverage::FileCoverage(const CoverageSpace & spac, const CodeFile & file, const TestSet & tset)
	: space(spac), cfile(file), tests(tset), lines() {}
void FileCoverage::load_in() {
	/* validation */
	if (cfile.get_text() == nullptr) {
		CError error(CErrorType::Runtime, "FileCoverage::load_in", "Code file is not loaded --> " + cfile.get_file().get_path());
		CErrorConsumer::consume(error);
	}
	else release();

	/* construct coverage for lines */
	const TextBuild & text = *(cfile.get_text());
	size_t line_number = text.numberOfLines(), k = 0;
	while (k++ < line_number) lines.push_back(nullptr);
}
void FileCoverage::release() {
	/* clear lines coverage */
	auto lbeg = lines.begin(), lend = lines.end();
	while (lbeg != lend) {
		LineCoverage * lcov = *(lbeg++);
		if (lcov != nullptr) delete lcov;
	}
	lines.clear();
}
void FileCoverage::cover(size_t line, TestCase::ID tid) {
	/* validation */
	if (line > lines.size() || line == 0) {
		CError error(CErrorType::Runtime, "FileCoverage::cover",
			"Invalid line (" + std::to_string(line)
			+ ") in file [" + std::to_string(lines.size()) + "]");
		CErrorConsumer::consume(error);
	}
	else if (!tests.has_test(tid)) return;

	/* get line-coverage */
	if (lines[line - 1] == nullptr) {
		lines[line - 1] = new LineCoverage(*this, line,
			tests.get_space().number_of_tests());
	}
	LineCoverage & lcov = *(lines[line - 1]);

	/* set the bit in line-coverage by tid */
	lcov.cover(tid);
}
bool FileCoverage::is_covered(const size_t line) const {
	if ((line != 0) && (line <= lines.size()))
		return lines[line - 1] != nullptr;
	else return false;
}
const LineCoverage & FileCoverage::get_line_coverage(const size_t line) const {
	/* validation */
	if (line > lines.size() || line == 0 || lines[line - 1] == nullptr) {
		CError error(CErrorType::Runtime, "FileCoverage::get_line_coverage",
			"Uncovered line (" + std::to_string(line)
			+ ") in file [" + std::to_string(lines.size()) + "]");
		CErrorConsumer::consume(error);
	}
	return *(lines[line - 1]);
}

CoverageSpace::CoverageSpace(const CTrace & p,
	const CoverageSource & src, const CodeSpace & cs, const TestSpace & ts)
	: project(p), source(src), code_space(cs), test_space(ts) {}
bool CoverageSpace::has_file_coverage(const CodeFile & cfile) const {
	const CodeSpace * cspace = &(cfile.get_space());
	if (cspace != &code_space) {
		CError error(CErrorType::InvalidArguments, "CoverageSpace::has_file_coverage",
			"Unmatched code space for file: \"" + cfile.get_file().get_path() + "\"");
		CErrorConsumer::consume(error);
	}
	return files.count(&cfile) > 0;
}
FileCoverage & CoverageSpace::get_file_coverage(const CodeFile & cfile) const {
	if (files.count(&cfile) == 0) {
		CError error(CErrorType::Runtime, "CoverageSpace::get_file_coverage",
			"Undefined code file: " + cfile.get_file().get_path());
		CErrorConsumer::consume(error);
	}
	auto iter = files.find(&cfile);
	return *(iter->second);
}
void CoverageSpace::add_file_coverage(const CodeFile & cfile, const TestSet & tests) {
	/* validation */
	const CodeSpace * cspace = &(cfile.get_space());
	const TestSpace * tspace = &(tests.get_space());
	if (cspace != &code_space) {
		CError error(CErrorType::InvalidArguments, "CoverageSpace::add_file_coverage",
			"Unmatched code space for file: \"" + cfile.get_file().get_path() + "\"");
		CErrorConsumer::consume(error);
	}
	else if (tspace != &test_space) {
		CError error(CErrorType::InvalidArguments, "CoverageSpace::add_file_coverage",
			"Unmatched test space for test: \"" + tspace->get_source().get_input_dir().get_path() + "\"");
		CErrorConsumer::consume(error);
	}

	/* remove original coverage */
	if (files.count(&cfile) > 0) {
		auto iter = files.find(&cfile);
		delete (iter->second);
		files.erase(&cfile);
	}
	files[&cfile] = new FileCoverage(*this, cfile, tests);
}
void CoverageSpace::del_file_coverage(const CodeFile & cfile) {
	if (files.count(&cfile) == 0) {
		CError error(CErrorType::InvalidArguments, "CoverageSpace::del_file_coverage",
			"Undefined code file: \"" + cfile.get_file().get_path() + "\"");
		CErrorConsumer::consume(error);
	}
	else {
		auto iter = files.find(&cfile);
		delete (iter->second);
		files.erase(&cfile);
	}
}
void CoverageSpace::clear() {
	auto beg = files.begin(), end = files.end();
	while (beg != end) {
		delete (beg++)->second;
	}
	files.clear();
}

void CoverageLoader::parse(FileCoverage & fcov) const {
	/* declarations */
	const TestSet & tests = fcov.get_tests();
	const CodeFile & cfile = fcov.get_file();
	const File & traces = fcov.get_space().get_source().get_trace_dir();
	const std::string & cfile_name = cfile.get_file().get_local_name();

	/* get ../traces/xxx.c/ */
	const File * cfile_trace = nullptr;
	auto files = traces.list_files();
	auto beg = files.begin(), end = files.end();
	while (beg != end) {
		File * file = *(beg++);

		if (file->is_directory()) {
			const std::string & name = file->get_local_name();
			if (name == cfile_name) {
				cfile_trace = file; break;
			}
		}
	}

	/* iterate all sub-files to extract coverage from ../traces/xxx.c/* (gcov) files */
	if (cfile_trace != nullptr) {
		auto tfiles = cfile_trace->list_files();
		auto tbeg = tfiles.begin(), tend = tfiles.end();
		while (tbeg != tend) {
			/* get the next gcov file */
			const File & tfile = *(*(tbeg++));
			if (tfile.is_directory()) continue;

			/* get the test case id */
			std::string name = tfile.get_local_name();
			if (name[0] != 't') continue;
			else name = name.substr(1, name.length() - 1);
			TestCase::ID tid = std::stoul(name);

			/* extract coverage from gcov file into FileCoverage */
			if (tests.has_test(tid))
				parse(fcov, tid, tfile);
		} /* end while: ../traces/xxx.c/* */
	}
}
void CoverageLoader::parse(FileCoverage & fcov, TestCase::ID tid, const File & tfile) const {
	/* declarations */
	LineReader reader(tfile.get_path());
	std::string linestr, numstr; char ch;
	int n, k; size_t line, count;

	/* iterate the lines in gcov file */
	while (reader.hasNext()) {
		/* get next line */
		linestr = reader.next();

		/* initialization */
		line = 0, count = 0;
		k = 0, n = linestr.length();

		/* derive the string for counts */
		{
			numstr = "";
			while (k < n) {
				ch = linestr[k++];
				if (ch == ':') break;
				else if (ch >= '0' && ch <= '9')
					numstr += ch;
			}
			if (numstr.empty()) continue;
			else count = std::stoul(numstr);
		}

		/* derive the string for line */
		{
			numstr = "";
			while (k < n) {
				ch = linestr[k++];
				if (ch == ':') break;
				else if (ch >= '0' && ch <= '9')
					numstr += ch;
			}
			if (numstr.empty()) continue;
			else line = std::stoul(numstr);
		}

		/* set coverage for each line */
		if (line > 0 && count > 0)
			fcov.cover(line, tid);
	} /* end while: iterate lines in xxx.gcov */
}

CTrace::CTrace(const File & root, const CodeSpace &cspace, const TestSpace &tspace) : dir(root), loader() {
	/* validation */
	if (!root.is_directory()) {
		CError error(CErrorType::InvalidArguments, "CTrace::CTrace", root.get_path() + " is not directory!");
		CErrorConsumer::consume(error);
	}

	/* find the ../traces/ */
	const File * traces = nullptr;
	auto files = root.list_files();
	auto beg = files.begin(), end = files.end();
	while (beg != end) {
		File * file = *(beg++);
		const std::string & name = file->get_local_name();
		if (name == "traces") {
			traces = file; break;
		}
	}

	/* create source */
	if (traces == nullptr) {
		CError error(CErrorType::Runtime, "CTrace::CTrace",
			"/traces is not found: \"" + root.get_path() + "\"");
		CErrorConsumer::consume(error);
	}
	else source = new CoverageSource(*traces);

	/* create space */
	space = new CoverageSpace(*this, *source, cspace, tspace);
}
bool CTrace::load_coverage(const CodeFile & cfile) {
	FileCoverage & fcov = space->get_file_coverage(cfile);
	loader.release(fcov);
	loader.load_in(fcov);
	loader.parse(fcov);
	return true;
}
