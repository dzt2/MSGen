#include "domset.h"

void ScoreMatrix::add_score_vectors(
	ScoreProducer & producer, ScoreConsumer & consumer) {
	ScoreVector * svec; 
	while ((svec = producer.produce()) != nullptr) {
		/* get the next score vector */
		Mutant::ID mid = svec->get_mutant();
		const BitSeq & bits = svec->get_vector();

		/* update the matrix */
		BitSeq::size_t bias = 
			mid * tspace.number_of_tests();
		BitSeq::size_t i, n = bits.bit_number();
		for (i = 0; i < n; i++) {
			bit bi = bits.get_bit(i);
			matrix->set_bit(bias + i, bi);
		}

		/* delete score vector */
		consumer.consume(svec); 
	}
}

void DomSetBuilder::derive_score_set(Mutant::ID mid, std::set<TestCase::ID> & scoreset) {
	TestSpace & tspace = matrix->get_test_space();
	TestCase::ID tid, n = tspace.number_of_tests();

	scoreset.clear();
	for (tid = 0; tid < n; tid++) {
		if (matrix->get_result(mid, tid))
			scoreset.insert(tid);
	}
}
bool DomSetBuilder::is_killed_by_all(Mutant::ID mj, const std::set<TestCase::ID> & scoreset) {
	compare_counts = compare_counts + 1;		// get the efficiency analysis

	auto beg = scoreset.begin();
	auto end = scoreset.end();
	while (beg != end) {
		TestCase::ID tid = *(beg++);
		if (!(matrix->get_result(mj, tid)))
			return false;
	}
	return true;
}

void DomSetBuilder_Greedy::erase_subsummeds(Mutant::ID mi, std::set<Mutant::ID> & M) {
	M.erase(mi);		// push

	/* get score set of mi */
	std::set<TestCase::ID> scoreset;
	derive_score_set(mi, scoreset);
	if (scoreset.empty()) return;

	/* compute those subsumed by mi */
	std::set<Mutant::ID> erases;
	auto beg = M.begin(), end = M.end();
	while (beg != end) {
		Mutant::ID mj = *(beg++);
		if (is_killed_by_all(mj, scoreset))
			erases.insert(mj);
	}

	/* update the set M */
	beg = erases.begin();
	end = erases.end();
	while (beg != end)
		M.erase(*(beg++));

	M.insert(mi);		// pope
}
bool DomSetBuilder_Greedy::get_next_mutants(
	const std::set<Mutant::ID> & M, 
	const std::set<Mutant::ID> & records,
	Mutant::ID & next) {
	auto beg = M.begin(), end = M.end();
	while (beg != end) {
		next = *(beg++);
		if (records.count(next) == 0)
			return true;
	}
	return false;
}
void DomSetBuilder_Greedy::compute(MutSet & ans) {
	/* declarations */
	std::set<Mutant::ID> records;
	std::set<Mutant::ID> domset;

	/* initialization */
	MutantSpace & mspace = matrix->get_mutant_space();
	Mutant::ID msize = mspace.number_of_mutants(), mid;
	for (mid = 0; mid < msize; mid++) 
		domset.insert(mid);

	/* compute the dominator set by elimination-greedly */
	Mutant::ID mi;
	while (get_next_mutants(domset, records, mi)) {
		records.insert(mi);
		erase_subsummeds(mi, domset);
	}

	/* update the answer */
	auto beg = domset.begin();
	auto end = domset.end();
	ans.clear_mutants();
	while (beg != end)
		ans.add_mutant(*(beg++));
}

