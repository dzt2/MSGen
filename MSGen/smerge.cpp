#include "smerge.h"

void SuMutantSet::update_mutants() {
	mutants->clear();	/* initialize the mutant set */

	/* get the nodes in graph that subsume all the others (may be equivalent) */
	const std::set<MuCluster *> & roots = graph.get_roots();
	auto beg = roots.begin(), end = roots.end();

	/* compute the equivalent cluster (if) */
	MuCluster * eqs = nullptr;
	if (roots.size() == 1) {
		MuCluster * cluster = *(beg++);
		if (cluster->get_score_degree() == 0)
			eqs = cluster;
	}
	
	/* nodes directly subsumed by equivalent cluster are subsuming */
	if (eqs != nullptr) {
		const std::vector<MuSubsume> & edges 
			= eqs->get_ou_port().get_edges();
		auto ebeg = edges.begin(), eend = edges.end();
		while (ebeg != eend) {
			const MuSubsume & edge = *(ebeg++);
			MuCluster & target = edge.get_target();
			const MutantSet & tmuts = target.get_mutants();
			mutants->disjunct(tmuts);
		}
	}
	/* otherwise, the roots in graph are subsuming mutants */
	else {
		beg = roots.begin(), end = roots.end();
		while (beg != end) {
			MuCluster * cluster = *(beg++);
			const MutantSet & tmuts = cluster->get_mutants();
			mutants->disjunct(tmuts);
		}
	}
}

void SuMutantMerger::open(SuMutantSet & result) {
	answer = &result;

	MSGraph & graph 
		= result.get_graph();
	builder.open(graph);
}
void SuMutantMerger::append(SuMutantSet & smuts) {
	/* get the subsuming mutants and its graph */
	const MutantSet & mutants = smuts.get_subsuming_mutants();
	Mutant::ID mid = 0, num = mutants.number_of_mutants();
	
	MSGraph & graph = smuts.get_graph();
	while (mid < num) {
		if (mutants.has_mutant(mid)) {
			MuCluster & cluster = graph.get_cluster_of(mid);
			const BitSeq & bits = cluster.get_score_vector();
			builder.add(mid, bits);
		}
		/* to the next mutant */ mid = mid + 1;
	}
}
void SuMutantMerger::extract() {
	builder.link(); answer->update_mutants();
}

const MutantSet & SuMutantExperimentCore::get_subsuming_mutants(const std::string & key) const {
	if (mut_map.count(key) == 0) {
		CError error(CErrorType::InvalidArguments, 
			"SuMutantExperimentCore::get_subsuming_mutants", 
			"Invalid key: \"" + key + "\"");
		CErrorConsumer::consume(error);
		exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = mut_map.find(key);
		SuMutantSet * smut = iter->second;
		return smut->get_subsuming_mutants();
	}
}
SuMutantSet & SuMutantExperimentCore::get_mutants(const std::string & key) {
	if (mut_map.count(key) == 0) {
		MSGraph * graph = new MSGraph(mspace);
		SuMutantSet * smut = new SuMutantSet(*graph);
		mut_map[key] = smut; keys.insert(key);
		return *smut;
	}
	else {
		auto iter = mut_map.find(key);
		SuMutantSet * smut = iter->second;
		return *smut;
	}
}
void SuMutantExperimentCore::clear_mutants() {
	auto beg = mut_map.begin();
	auto end = mut_map.end();
	while (beg != end) {
		SuMutantSet * smut = (beg++)->second;
		MSGraph & graph = smut->get_graph();
		delete smut; delete &graph;
	}
	mut_map.clear(); keys.clear();
}

void SuMutantExperimentDriver::derive_operator_I(
	ScoreProducer & producer, ScoreConsumer & consumer) {
	/* builders for initializing MSG according to operator */
	std::map<std::string, MSGBuilder *> builders;
	MutantSpace & mspace = core->get_space();

	/* create ulinker graph for mutants of each operator */
	ScoreVector * score_vector;
	while ((score_vector = producer.produce()) != nullptr) {
		/* get the next mutant */
		const BitSeq & bits = score_vector->get_vector();
		Mutant::ID mid = score_vector->get_mutant();
		Mutant & mutant = mspace.get_mutant(mid);
		const std::string & oprt = mutant.get_operator();

		/* get the subsuming set for the operator */
		SuMutantSet & smut = core->get_mutants(oprt);
		MSGraph & graph = smut.get_graph();

		/* get the builder */
		MSGBuilder * builder;
		if (builders.count(oprt) == 0) {
			builder = new MSGBuilder();
			builder->open(graph);
			builders[oprt] = builder;
		}
		else {
			auto iter = builders.find(oprt);
			builder = iter->second;
		}

		/* append the mutant into the graph-builder */
		builder->add(mid, bits);

		/* consume the vector */
		consumer.consume(score_vector);
	} /* end while */

	/* construct linked graph */
	auto beg = builders.begin(), end = builders.end();
	while(beg != end) {
		/* get the next operator and its subsuming mutants */
		const std::string & oprt = beg->first;
		MSGBuilder * builder = beg->second;
		SuMutantSet & smut = core->get_mutants(oprt);

		/* link the graph and delete builder */
		builder->link(); 
		builder->close(); 
		delete builder;

		/* update the subsuming mutants */
		smut.update_mutants();

		beg++;	/* to the next operator */
	}
	builders.clear();
}
void SuMutantExperimentDriver::derive_operator_II() {
	CError error(CErrorType::Runtime, 
		"SuMutantExperimentDriver::derive_operator_II", 
		"Invalid access: operators have not been designed");
	CErrorConsumer::consume(error);
	exit(CErrorType::Runtime);
}
void SuMutantExperimentDriver::derive_global_III() {
	CError error(CErrorType::Runtime,
		"SuMutantExperimentDriver::derive_global_III",
		"Invalid access: operators have not been designed");
	CErrorConsumer::consume(error);
	exit(CErrorType::Runtime);
}