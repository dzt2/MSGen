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

void DomSetBuilder_Greedy::derive_score_set(Mutant::ID mid, std::set<TestCase::ID> & scoreset) {
	TestSpace & tspace = matrix->get_test_space();
	TestCase::ID tid, n = tspace.number_of_tests();

	scoreset.clear();
	for (tid = 0; tid < n; tid++) {
		if (matrix->get_result(mid, tid))
			scoreset.insert(tid);
	}
}
bool DomSetBuilder_Greedy::is_killed_by_all(Mutant::ID mj, const std::set<TestCase::ID> & scoreset) {
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
	MutantSpace & mspace = ans.get_space();
	Mutant::ID msize = mspace.number_of_mutants(), mid;
	for (mid = 0; mid < msize; mid++) domset.insert(mid);

	/* compute the dominator set by elimination-greedly */
	Mutant::ID mi; compare_counts = 0;
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
	compare_counts++;

	auto beg = tests.begin();
	auto end = tests.end();
	while (beg != end) {
		TestCase::ID tid = *(beg++);
		if (!(matrix->get_result(mid, tid)))
			return false;
	}
	return true;
}
void DomSetBuilder_Blocks::erase_subsummeds(Mutant::ID mid,
	std::set<Mutant::ID> & M, const std::set<TestCase::ID> & scoreset) {
	M.erase(mid);						// push

	if (!scoreset.empty()) {
		/* mutants to be eliminated */
		std::set<Mutant::ID> trash;

		/* compute those to be subsumed by mid */
		auto beg = M.begin();
		auto end = M.end();
		while (beg != end) {
			Mutant::ID mid = *(beg++);
			if (is_killed_by_all(mid, scoreset))
				trash.insert(mid);
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
		records.insert(mid);	// record the mutant as visited

		/* eliminate those subsummed by mid */
		const std::set<Mutant::ID> & scoreset 
			= derive_score_set(coverset, mid);
		erase_subsummeds(mid, domset, scoreset);
	} // end while
	records.clear();
}

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
			MSG_Node & block = graph.get_node_of(mid);
			if (blocks.count(&block)) mutants.insert(mid);
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
	while (beg != end) {
		Mutant::ID mi = *(beg++);

		/* get feasible domain of the mutant */
		const std::set<TestCase::ID> & scoreset = derive_score_set(coverset, mi);
		std::set<MSG_Node *> domain; get_feasible_domain(scoreset, domain);
		domain.erase(&block);

		/* get the mutants in feasible domain from M */
		std::set<Mutant::ID> mutants; collect_mutants_of(M, domain, mutants);

		/* eliminate those subsumed by mi in mutants */
		erase_subsummeds(mi, mutants, scoreset);

		/* reput the remainders into the set M */
		auto beg2 = mutants.begin(), end2 = mutants.end();
		while (beg2 != end2) M.insert(*(beg2++));
	} // emd wjo;e
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
	MSG_Node * root = nullptr;
	MutantSpace & mspace = ans.get_space();
	Mutant::ID msize = mspace.number_of_mutants();
	for (Mutant::ID mid = 0; mid < msize; mid++) {
		if (graph.has_node_of(mid)) {
			MSG_Node & block = graph.get_node_of(mid);
			if (block.get_score_degree() > 0) 
				M.insert(mid);
			if (block.get_in_port().degree() == 0)
				root = &block;
		}
	}

	/* compute by the block from top to down */
	std::queue<MSG_Node *> tqueue; 
	std::queue<MSG_Node *> list; 
	std::set<MSG_Node *> records;
	tqueue.push(root);
	while (!tqueue.empty()) {
		MSG_Node * block = 
			tqueue.front(); tqueue.pop(); 

		if(block->get_score_degree() > 0)
			list.push(block);

		const MSG_Port & port = block->get_ou_port();
		for (int k = 0; k < port.degree(); k++) {
			MSG_Edge & edge = port.get_edge(k);
			MSG_Node & trg = edge.get_target();
			if (records.count(&trg) == 0) {
				records.insert(&trg);
				tqueue.push(&trg);
			}
		}
	}

	/* compute by each block */
	while (!list.empty()) {
		MSG_Node & block = *(list.front());
		update_by_block(M, block); list.pop();
	}
}

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






