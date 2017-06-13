#include "mclass.h"

MuClass::MuClass(MuClassSet & cset, MuFeature ft) : classes(cset), feature(ft) {
	mutants = classes.get_mutants().get_space().create_set();
}
MuClass::~MuClass() { classes.get_mutants().get_space().delete_set(mutants); }

MuClass & MuClassSet::get_class(MuFeature ft) const {
	if (classes.count(ft) == 0) {
		CError error(CErrorType::InvalidArguments, "MuClassSet::get_class", "Undefined class");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = classes.find(ft); return *(iter->second);
	}
}
MuClassSet::~MuClassSet() { clear_classes(); }
MuClass * MuClassSet::new_class(MuFeature ft) {
	if (classes.count(ft) > 0) {
		CError error(CErrorType::InvalidArguments, "MuClassSet::new_class", "Duplicated feature");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		MuClass * _class = new MuClass(*this, ft);
		classes[ft] = _class; return _class;
	}
}
void MuClassSet::clear_classes() {
	auto beg = classes.begin();
	auto end = classes.end();
	while (beg != end)
		delete (beg++)->second;
	classes.clear();
}

void MuClassifier::classify(MuClassSet & class_set) {
	/* to classify mutants */
	Mutant::ID mid = 0; MuFeature ft; class_set.clear_classes();
	const MutantSet & mutants = class_set.get_mutants();
	size_t num = mutants.get_space().number_of_mutants();
	while (true) {
		/* get the next feature */ next(mutants, mid, ft);

		/* validation */
		if (mid >= num) break;
		else if (ft == nullptr) {
			CError error(CErrorType::Runtime, "MuClassifier::classify", "Classify Error");
			CErrorConsumer::consume(error); exit(CErrorType::Runtime);
		}
		
		/* get the class */
		MuClass * _class;
		if (class_set.has_class(ft))
			 _class = &(class_set.get_class(ft));
		else _class = class_set.new_class(ft);

		/* insert mutant to the class */
		_class->add_mutant(mid);

		/* roll to the next mutant */ mid++;
	} /* end while */

	/* return */ return;
}
MuClassifier::~MuClassifier() {}

void MuClassifierByOperator::next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & feature) {
	size_t num = mutants.get_space().number_of_mutants();
	while (mid < num && !mutants.has_mutant(mid)) mid++;

	if (mid < num) {
		/* get the mutant operator */
		Mutant & mutant = mutants.get_space().get_mutant(mid);
		const std::string & oprt = mutant.get_operator();

		/* get the mutant feature */
		if (operators.count(oprt) == 0) {
			feature = new std::string(oprt);
			operators[oprt] = (std::string *) feature;
		}
		else {
			auto iter = operators.find(oprt);
			feature = iter->second;
		}
	}
	else feature = nullptr;
}
void MuClassifierByLocation::next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & feature) {
	/* get next feature */
	size_t num = mutants.get_space().number_of_mutants();
	while (mid < num && !mutants.has_mutant(mid)) mid++;

	if (mid < num) {
		/* get next mutant and its location (text) */
		Mutant & mutant = mutants.get_space().get_mutant(mid);
		const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & loc = mutation.get_location();
		std::string loctext = std::to_string(loc.get_bias());
		loctext += ":" + std::to_string(loc.get_length());

		/* get the feature */
		if (locations.count(loctext) == 0) {
			feature = new CodeLocation(loc);
			locations[loctext] = (CodeLocation *)feature;
		}
		else {
			auto iter = locations.find(loctext);
			feature = iter->second;
		}
	}
	else feature = nullptr;
}
void MuClassifierByCoverage::next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & feature) {
	/* initialization */
	mid = mutants.get_space().number_of_mutants();
	CoverageVector * covvec; feature = nullptr;
	
	/* iterate all coverage in the producer */
	while ((covvec = producer->produce()) != nullptr) {
		/* get the next mutant's feature */
		if (mutants.has_mutant(covvec->get_mutant())) {
			/* get the mutant id */ mid = covvec->get_mutant();

			/* get the leaf for this coverage */
			BitTrie * leaf = trie->insert_vector(covvec->get_coverage());

			/* first time to create coverage vector */
			if (leaf->get_data() == nullptr) {
				feature = new BitSeq(covvec->get_coverage());
				leaf->set_data(feature);
			}
			else feature = leaf->get_data();
		}

		/* consume the coverage */
		consumer->consume(covvec);
		/* has-mutant, get-feature, break */
		if (feature != nullptr) break;
	} /* end while */
}
void MuClassifierByScore::next(const MutantSet & mutants, Mutant::ID & mid, MuFeature & feature) {
	/* initialization */
	mid = mutants.get_space().number_of_mutants();
	ScoreVector * vec; feature = nullptr;

	/* iterate all coverage in the producer */
	while ((vec = producer->produce()) != nullptr) {
		/* get the next mutant's feature */
		if (mutants.has_mutant(vec->get_mutant())) {
			/* get the mutant id */ mid = vec->get_mutant();

			/* get the leaf for this coverage */
			BitTrie * leaf = trie->insert_vector(vec->get_vector());

			/* first time to create coverage vector */
			if (leaf->get_data() == nullptr) {
				feature = new BitSeq(vec->get_vector());
				leaf->set_data(feature);
			}
			else feature = leaf->get_data();
		}

		/* consume the coverage */
		consumer->consume(vec);
		/* has-mutant, get-feature, break */
		if (feature != nullptr) break;
	} /* end while */
}
