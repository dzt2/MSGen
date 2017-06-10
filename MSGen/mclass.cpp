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
MuClassSet::~MuClassSet() {
	auto beg = classes.begin();
	auto end = classes.end();
	while (beg != end)
		delete (beg++)->second;
	classes.clear();
}
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

MuClassSet & MuClassifier::classify(const MutantSet & mutants) {
	/* create a new class-set */
	MuClassSet * class_set = new MuClassSet(mutants);
	this->pool.insert(class_set);

	/* to classify mutants */
	Mutant::ID mid = 0; MuFeature ft;
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
		if (class_set->has_class(ft))
			 _class = &(class_set->get_class(ft));
		else _class = class_set->new_class(ft);

		/* insert mutant to the class */
		_class->add_mutant(mid);

		/* roll to the next mutant */ mid++;
	} /* end while */

	/* return */ return *class_set;
}
void MuClassifier::delete_classes(MuClassSet & _classes) {
	if (pool.count(&_classes) == 0) {
		CError error(CErrorType::InvalidArguments, "MuClassifier::delete_classes", "Undefined class-set");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		pool.erase(&_classes); delete &_classes;
	}
}
MuClassifier::~MuClassifier() {
	auto beg = pool.begin();
	auto end = pool.end();
	while (beg != end)
		delete *(beg++);
	pool.clear();
}

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
			BitTrie * leaf = trie.get_leaf(covvec->get_coverage());

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
			BitTrie * leaf = trie.get_leaf(vec->get_vector());

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

void test_classify_operator(const MutantSet & mutants) {
	// classifier
	MuClassifierByOperator classifier;
	MuClassSet & class_set = classifier.classify(mutants);

	// 
	const std::map<MuFeature, MuClass *>
		classes = class_set.get_classes();
	std::cout << "There are " << classes.size() << " classes generated...\n";

	auto beg = classes.begin(), end = classes.end();
	while (beg != end) {
		/* get next feature and class */
		MuFeature feature = beg->first;
		MuClass & _class = *(beg->second);
		beg++;

		/* for operator */
		std::string * oprt = (std::string *) feature;
		std::cout << "\t" << *oprt << " : " << _class.size() << "\n";
	}
}
void test_classify_location(const MutantSet & mutants) {
	MuClassifierByLocation classifier;
	MuClassSet & class_set = classifier.classify(mutants);

	const std::map<MuFeature, MuClass *> & classes = class_set.get_classes();
	std::cout << "There are " << classes.size() << " location-classes\n";
	
	auto beg = classes.begin(), end = classes.end();
	while (beg != end) {
		/* get next class and its location */
		MuClass & _class = *((beg++)->second);
		MuFeature feature = _class.get_feature();
		CodeLocation * loc = (CodeLocation *)feature;

		/* print location */
		std::cout << "\t[" << loc->get_bias() << ", " << loc->get_length() << "]\t" << _class.size() << "\n";
	}
}

int main() {
	// initialization
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string prname = "mid";
	File & root = *(new File(prefix + prname)); TestType ttype = TestType::general;

	// create code-project, mutant-project, test-project
	CProgram & program = *(new CProgram(root));
	CMutant & cmutant = *(new CMutant(root, program.get_source()));

	// load mutations 
	const CodeSpace & cspace = cmutant.get_code_space();
	const std::set<CodeFile *> & cfiles = cspace.get_code_set();
	auto cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		const CodeFile & cfile = *(*(cfile_beg++));
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		cmutant.load_mutants_for(mspace, true);
		std::cout << "Load " << mspace.number_of_mutants() <<
			" mutants for: " << cfile.get_file().get_path() << "\n" << std::endl;
	}

	// classify
	cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		// get next file and load its text 
		CodeFile & cfile = *(*(cfile_beg++)); cspace.load(cfile);

		// get mutations for code-file
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);
		MutantSet & mutants = *(mspace.create_set());
		mutants.complement();

		// test 
		test_classify_location(mutants);
	}

	// delete resources
	delete &cmutant; delete &program; delete &  root;

	// exit 
	std::cout << "\nPress any key to exit..."; getchar(); exit(0);
}