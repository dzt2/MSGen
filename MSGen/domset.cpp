#include "domset.h"
#include <cstdlib>
#include <time.h>

/// static tool method
static void random_sort(std::vector<Mutant::ID> & list) {
	size_t limit = list.size();
	while (limit > 1) {
		/* random index */
		srand((unsigned)time(nullptr));
		int next = rand() % limit;
		/* swap */
		auto temp = list[next];
		list[next] = list[limit - 1];
		list[limit - 1] = temp;
		/* to next */ limit--;
	}
}
static void random_sort(std::vector<MSG_Node *> & list) {
	size_t limit = list.size();
	while (limit > 1) {
		/* random index */
		srand((unsigned)time(nullptr));
		int next = rand() % limit;
		/* swap */
		auto temp = list[next];
		list[next] = list[limit - 1];
		list[limit - 1] = temp;
		/* to next */ limit--;
	}
}

/// score matrix
void ScoreMatrix::add_score_vectors(
	ScoreProducer & producer, ScoreConsumer & consumer) {
	ScoreVector * svec; 
	equivalents = 0;
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

		if (svec->get_degree() == 0)
			equivalents++;

		/* delete score vector */
		consumer.consume(svec); 
	}
}

/// greedy algorithm
void DomSetBuilder_Greedy::derive_score_set(Mutant::ID mid, std::set<TestCase::ID> & scoreset) {
	TestSpace & tspace = matrix->get_test_space();
	TestCase::ID tid, n = tspace.number_of_tests();

	scoreset.clear();
	for (tid = 0; tid < n; tid++) {
		if (matrix->get_result(mid, tid))
			scoreset.insert(tid);
		counter->testing(mid, tid);
	}
}
bool DomSetBuilder_Greedy::is_killed_by_all(Mutant::ID mj, const std::set<TestCase::ID> & scoreset) {
	auto beg = scoreset.begin();
	auto end = scoreset.end();
	while (beg != end) {
		TestCase::ID tid = *(beg++);
		counter->testing(mj, tid);
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

	/* determine equivalent */
	counter->equivalent(mi);
	if (scoreset.empty()) return;

	/* compute those subsumed by mi */
	std::set<Mutant::ID> erases;
	auto beg = M.begin(), end = M.end();
	while (beg != end) {
		Mutant::ID mj = *(beg++);
		counter->compare(mi, mj);
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
	auto beg = M.begin();
	auto end = M.end();
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
	MutantSpace & mspace = ans.get_space();
	Mutant::ID msize = mspace.number_of_mutants(), mid;
	for (mid = 0; mid < msize; mid++) domset.insert(mid);

	/* compute the dominator set by elimination-greedly */
	Mutant::ID mi; 
	while (get_next_mutants(domset, records, mi)) {
		records.insert(mi);
		counter->decline(mi, domset.size());
		erase_subsummeds(mi, domset);
	}

	/* update the answer */
	auto beg = domset.begin();
	auto end = domset.end();
	ans.clear_mutants();
	while (beg != end)
		ans.add_mutant(*(beg++));
}

/// coverage-based algorithm
/// inner-block computations
void DomSetBuilder_Blocks::collect_mutant_of(const std::set<Mutant::ID> & M, MSG_Node & node, std::set<Mutant::ID> & mutants) {
	mutants.clear();

	auto beg = M.begin();
	auto end = M.end();
	while (beg != end) {
		Mutant::ID mid = *(beg++);
		if (node.has_mutant(mid))
			mutants.insert(mid);
	}
}
const std::set<TestCase::ID> & DomSetBuilder_Blocks::derive_cover_set(MSG_Node & block) {
	/* first time */
	if (cover_sets.count(&block) == 0) {
		/* create a new coverage set for the block in cache */
		std::set<TestCase::ID> & coverset 
			= *(new std::set<TestCase::ID>());
		cover_sets[&block] = &coverset;

		/* extract test id from the block coverage bits */
		const BitSeq & bits = block.get_score_vector();
		TestSpace & tspace = matrix->get_test_space();
		TestCase::ID tsize = tspace.number_of_tests();
		for (TestCase::ID tid = 0; tid < tsize; tid++)
			if (bits.get_bit(tid)) coverset.insert(tid);

		return coverset;	// return;
	}
	/* has been computed */
	else {
		auto iter = cover_sets.find(&block);
		return *(iter->second);
	}
}
const std::set<TestCase::ID> &  DomSetBuilder_Blocks::derive_score_set(
	const std::set<TestCase::ID> & coverset, Mutant::ID mid) {
	/* first time */
	if (score_sets.count(mid) == 0) {
		/* create a new score set for the mutant in cache */
		std::set<TestCase::ID> & scoreset
			= *(new std::set<TestCase::ID>());
		score_sets[mid] = &scoreset;

		/* get the tests that kill mutant within coverage-set */
		auto beg = coverset.begin();
		auto end = coverset.end();
		while (beg != end) {
			TestCase::ID tid = *(beg++);
			counter->testing(mid, tid);
			if (matrix->get_result(mid, tid))
				scoreset.insert(tid);
		}

		return scoreset;	// return 
	}
	/* has been computed */
	else {
		auto iter = score_sets.find(mid);
		return *(iter->second);
	}
}
bool DomSetBuilder_Blocks::is_killed_by_all(Mutant::ID mid, const std::set<TestCase::ID> & tests) {
	auto beg = tests.begin();
	auto end = tests.end();
	while (beg != end) {
		TestCase::ID tid = *(beg++);
		counter->testing(mid, tid);
		if (!(matrix->get_result(mid, tid)))
			return false;
	}
	return true;
}
void DomSetBuilder_Blocks::erase_subsummeds(Mutant::ID mid,
	std::set<Mutant::ID> & M, const std::set<TestCase::ID> & scoreset) {
	M.erase(mid);						// push

	counter->equivalent(mid);
	if (!scoreset.empty()) {
		/* mutants to be eliminated */
		std::set<Mutant::ID> trash;

		/* compute those to be subsumed by mid */
		auto beg = M.begin();
		auto end = M.end();
		while (beg != end) {
			Mutant::ID mj = *(beg++);
			counter->compare(mid, mj);
			if (is_killed_by_all(mj, scoreset))
				trash.insert(mj);
		}

		/* eliminate those subsummed by mid */
		beg = trash.begin();
		end = trash.end();
		while (beg != end) 
			M.erase(*(beg++));
	}
	else return;		// equivalent mutants are not used here

	M.insert(mid);						// pope
}
bool DomSetBuilder_Blocks::get_next_mutants(const std::set<Mutant::ID> & M,
	const std::set<Mutant::ID> & records, Mutant::ID & next) {
	auto beg = M.begin();
	auto end = M.end();
	while (beg != end) {
		next = *(beg++);
		if (records.count(next) == 0)
			return true;
	}
	return false;
}
void DomSetBuilder_Blocks::compute_inner_block(std::set<Mutant::ID> & M, 
	MSG_Node & block, std::set<Mutant::ID> & domset) {
	/* initialization */
	collect_mutant_of(M, block, domset);
	auto beg = domset.begin();
	auto end = domset.end();
	while (beg != end) M.erase(*(beg++));

	/* collect the coverage set */
	const std::set<TestCase::ID> & 
		coverset = derive_cover_set(block);
	
	/* erase redundant mutants */
	Mutant::ID mid;
	std::set<Mutant::ID> records;
	while (get_next_mutants(domset, records, mid)) {
		// record the mutant as visited
		records.insert(mid);	
		counter->decline(mid, domset.size());

		/* get the score set of this mutant */
		const std::set<Mutant::ID> & scoreset
			= derive_score_set(coverset, mid);

		/* eliminate those subsummed by mid */
		erase_subsummeds(mid, domset, scoreset);
	} // end while
	records.clear();
}
/// inter-block computations
void DomSetBuilder_Blocks::get_feasible_domain(const 
	std::set<TestCase::ID> & scoreset, std::set<MSG_Node *> & domain) {
	domain.clear();				// initialization

	size_t bsize = graph.size();
	for (size_t k = 0; k < bsize; k++) {
		/* get the next block */
		MSG_Node & block = graph.get_node(k);
		const BitSeq & bits = block.get_score_vector();
		
		/* whether the block is feasible */
		bool is_feasible = true;
		auto beg = scoreset.begin();
		auto end = scoreset.end();
		while (beg != end) {
			TestCase::ID tid = *(beg++);
			if (!bits.get_bit(tid)) {
				is_feasible = false;
				break;
			}
		}

		/* insert to the domain */
		if (is_feasible)
			domain.insert(&block);
	}
}
void DomSetBuilder_Blocks::collect_mutants_of(const std::set<Mutant::ID> & M,
	const std::set<MSG_Node *> & blocks, std::set<Mutant::ID> & mutants) {
	mutants.clear();

	auto beg = M.begin();
	auto end = M.end();
	while (beg != end) {
		Mutant::ID mid = *(beg++);
		if (graph.has_node_of(mid)) {
			MSG_Node & block = 
				graph.get_node_of(mid);
			if (blocks.count(&block) > 0) 
				mutants.insert(mid);
		}
	}
}
void DomSetBuilder_Blocks::compute_inter_block(std::set<Mutant::ID> & M, 
	MSG_Node & block, const std::set<Mutant::ID> domset) {
	/* get the coverage set as basis */
	const std::set<TestCase::ID> & 
		coverset = derive_cover_set(block);

	auto beg = domset.begin();
	auto end = domset.end();
	int bias = 0;
	while (beg != end) {
		Mutant::ID mi = *(beg++);
		counter->decline(mi, M.size() + domset.size() - (bias--));

		/* get feasible domain of the mutant */
		const std::set<TestCase::ID> & scoreset = derive_score_set(coverset, mi);
		std::set<MSG_Node *> domain; get_feasible_domain(scoreset, domain);
		domain.erase(&block);

		/* get the mutants in feasible domain from M */
		std::set<Mutant::ID> mutants; collect_mutants_of(M, domain, mutants);
		auto beg1 = mutants.begin(), end1 = mutants.end();
		while (beg1 != end1) M.erase(*(beg1++));

		/* eliminate those subsumed by mi in mutants */
		erase_subsummeds(mi, mutants, scoreset);

		/* reput the remainders into the set M */
		auto beg2 = mutants.begin(), end2 = mutants.end();
		while (beg2 != end2) M.insert(*(beg2++));
	} // emd wjo;e
}
/// coverage-based implementation
MSG_Node * DomSetBuilder_Blocks::get_min_mutants(const std::set<MSG_Node *> & blocks) {
	size_t number = UINT32_MAX;
	auto beg = blocks.begin();
	auto end = blocks.end();

	MSG_Node * min = nullptr;
	while(beg != end) {
		MSG_Node * node = *(beg++);
		if (node->get_mutants().number_of_mutants() < number) {
			number = node->get_mutants().number_of_mutants();
			min = node;
		}
	}
	return min;
}
void DomSetBuilder_Blocks::sort_blocks_by(MS_Graph & cgraph, std::vector<MSG_Node *> & list) {
	list.clear();

	/* collect by degree */
	std::map<size_t, std::set<MSG_Node *> *> degree_blocks;
	std::vector<size_t> block_degrees; size_t n = cgraph.size();
	for (size_t k = 0; k < n; k++) {
		MSG_Node & block = cgraph.get_node(k);
		if (block.get_score_degree() > 0) {
			/* get the degree of this block */
			size_t degree = block.get_score_degree();

			/* insert for first time */
			if (degree_blocks.count(degree) == 0) {
				degree_blocks[degree] = new std::set<MSG_Node *>();
				block_degrees.push_back(degree);
			}

			/* update the block set for degree */
			auto iter = degree_blocks.find(degree);
			(iter->second)->insert(&block);
		}
	}

	/* sort by degree */
	std::sort(block_degrees.begin(), block_degrees.end());
	for (int i = 0; i < block_degrees.size(); i++) {
		/* get nodes for the degree */
		size_t degree = block_degrees[i];
		auto iter = degree_blocks.find(degree);
		std::set<MSG_Node *> & blocks = *(iter->second);

		/* get the list of sorted line */
		std::vector<MSG_Node *> vecs;
		auto beg = blocks.begin();
		auto end = blocks.end();
		while (beg != end)
			vecs.push_back(*(beg++));
		delete &blocks;

		/* randomly sort the nodes into list */
		random_sort(vecs);
		for (int k = 0; k < vecs.size(); k++)
			list.push_back(vecs[k]);
	}
	degree_blocks.clear(); block_degrees.clear();

}
void DomSetBuilder_Blocks::update_by_block(std::set<Mutant::ID> & M, MSG_Node & block) {
	if (block.get_score_degree() == 0) return;
	else {
		/* elimination */
		std::set<Mutant::ID> domset;
		compute_inner_block(M, block, domset);
		compute_inter_block(M, block, domset);

		/* reput dominator one */
		auto beg = domset.begin(), end = domset.end();
		while (beg != end) M.insert(*(beg++)); domset.clear();
	}
}
void DomSetBuilder_Blocks::compute(MutSet & ans) {
	/* initialization */
	std::set<Mutant::ID> M; 
	MutantSpace & mspace = ans.get_space();
	Mutant::ID msize = mspace.number_of_mutants();
	for (Mutant::ID mid = 0; mid < msize; mid++) {
		if (graph.has_node_of(mid)) {
			MSG_Node & block = graph.get_node_of(mid);
			if (block.get_score_degree() > 0) 
				M.insert(mid);
		}
	}

	/* sort block based on its dominance */
	std::vector<MSG_Node *> list;
	sort_blocks_by(graph, list);

	/* compute by each block */
	int n = list.size();
	for (int i = 0; i < n; i++) {
		MSG_Node & block = *(list[i]);
		update_by_block(M, block); 
	}

	/* put into dominator set */
	auto beg = M.begin();
	auto end = M.end();
	while (beg != end)
		ans.add_mutant(*(beg++));
}
/// resource methods
void DomSetBuilder_Blocks::clear_cache() {
	auto beg1 = cover_sets.begin();
	auto end1 = cover_sets.end();
	while (beg1 != end1) 
		delete ((beg1++)->second);

	auto beg2 = score_sets.begin();
	auto end2 = score_sets.end();
	while (beg2 != end2) 
		delete ((beg2++)->second);
}
