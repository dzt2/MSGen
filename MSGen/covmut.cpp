#include "cscore.h"
#include "ctrace.h"

/*
int main() {
	// declarations
	std::string name = "tcas";
	File & root = *(new File("../../SiemensSuite/" + name));
	CProgram & prog = *(new CProgram(root));
	CTest & tesp = *(new CTest(TestType::general, root, prog.get_exec()));
	CMutant & cmut = *(new CMutant(root, prog.get_source()));

	// load tests
	tesp.load(); const TestSpace &tspace = tesp.get_space();
	std::cout << "Load " << tspace.number_of_tests() << " tests...\n";
	// load mutants 
	const CodeSpace & cspace = prog.get_source();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto cbeg = cfiles.begin(), cend = cfiles.end();
	while (cbeg != cend) {
		const CodeFile & cfile = *(*(cbeg++));
		MutantSpace & mspace = cmut.get_mutants_of(cfile);
		cmut.load_mutants_for(mspace, true);
		std::cout << "Load " << mspace.number_of_mutants() <<
			" mutants at \"" << cfile.get_file().get_path() << "\"\n";
	}

	// get trace and score 
	CScore & score = *(new CScore(root, cmut, tesp));
	CTrace & trac = *(new CTrace(root, cspace, tspace));
	TestSet & tests = *(tesp.malloc_test_set()); tests.complement();

	// get score and coverage vector
	const CoverageSpace & covspac = trac.get_space();
	cbeg = cfiles.begin(), cend = cfiles.end();
	while (cbeg != cend) {
		// gets
		CodeFile & cfile = *(*(cbeg++));
		MutantSpace & mspace = cmut.get_mutants_of(cfile);
		MutantSet & mutants = *(mspace.create_set()); mutants.complement();
		CoverageSpace & covspac = trac.get_space();

		// load code text
		cspace.load(cfile);

		// load coverage
		covspac.add_file_coverage(cfile, tests);
		trac.load_coverage(cfile);
		std::cout << "Coverage file...\n";

		// get file coverage and coverage producer
		FileCoverage & filecov = covspac.get_file_coverage(cfile);
		CoverageProducer cproducer(mspace, filecov);
		CoverageConsumer cconsumer;

		// get score producer and coverage producer
		ScoreSource & score_src = score.get_source(cfile);
		ScoreFunction & func = *(score_src.create_function(tests, mutants));
		ScoreProducer producer(func); ScoreConsumer consumer(func);

		// reading
		int not_subsume = 0, strict_subsume = 0, equal = 0;
		CoverageVector * cvec; ScoreVector * svec;
		while ((cvec = cproducer.produce()) != nullptr) {
			svec = producer.produce();
			// if (svec == nullptr) break;

			// not synchronized
			if (svec->get_mutant() != cvec->get_mutant()) {
				std::cerr << "Error: " << svec->get_mutant() << " ~~~ " << cvec->get_mutant() << "\n";
				getchar();
			}

			if (svec != nullptr) {
				// get mutant 
				const BitSeq & cbits = cvec->get_coverage();
				const BitSeq & sbits = svec->get_vector();
				Mutant & mutant = mspace.get_mutant(cvec->get_mutant());
				const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
				const TextBuild & text = *(cfile.get_text());
				size_t line = text.lineOfIndex(mutation.get_location().get_bias());

				// some tests that kill mutant does not cover it
				if (!sbits.subsume(cbits)) {
					std::cerr << "Impossible for mutant #" << cvec->get_mutant() << "\n";
					std::cerr << "\tLine: " << line << "\n";
					std::cerr << "\tCover: " << cbits.to_string() << "\n";
					std::cerr << "\tKills: " << sbits.to_string() << "\n\n";
					not_subsume++;
				}
				// coverage vector strictly subsume score vector
				else if (!sbits.equals(cbits)) strict_subsume++;
				// equal with each other
				else equal++;
			}

			// delete both vectors 
			if (svec != nullptr)
				consumer.consume(svec);
			cconsumer.consume(cvec);
		} // end while

		std::cout << "\n[Incorrect: " << not_subsume << "; Strict: " << strict_subsume << "; Equals: " << equal << "]\n";

		// delete resources
		mspace.delete_set(&mutants);
		score_src.delete_function(&func);
	}

	// release resources
	delete &score; delete &trac; delete &tesp;
	delete &cmut; delete &prog; delete &root;
	std::cout << "Press any key to exit..." << std::endl;
	getchar(); return 0;
}
*/