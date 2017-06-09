#include "mclass.h"

const std::string & MutFeature::get_string() const {
	if (type != MutMetaType::String) {
		CError error(CErrorType::InvalidArguments, "MutFeature::get_string",
			"Invalid feature type: (" + std::to_string(type) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *((std::string *) content);
}
const CodeLocation & MutFeature::get_location() const {
	if (type != MutMetaType::Location) {
		CError error(CErrorType::InvalidArguments, "MutFeature::get_location",
			"Invalid feature type: (" + std::to_string(type) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *((CodeLocation *) content);
}
const BitSeq & MutFeature::get_bit_string() const {
	if (type != MutMetaType::BitSequence) {
		CError error(CErrorType::InvalidArguments, "MutFeature::get_bit_string",
			"Invalid feature type: (" + std::to_string(type) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else return *((BitSeq *)content);
}

MutClass::MutClass(MutClassGroup & grp, MutFeature & ft) 
	: group(grp), feature(&ft) {
	mutants = grp.get_mutant_space().create_set();
}
MutClass::~MutClass() {
	delete feature;
	group.get_mutant_space().delete_set(mutants);
}

MutClassGroup::MutClassGroup(MutantSpace & space)
	: mspace(space), classes(), index() {}
MutClass & MutClassGroup::get_class_for(Mutant::ID mid) const {
	if (index.count(mid) == 0) {
		CError error(CErrorType::InvalidArguments, "MutClassGroup::get_class_for",
			"Undefined mutant (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error); exit(0);
	}
	else {
		auto iter = index.find(mid); return *(iter->second);
	}
}
MutClass * MutClassGroup::new_class(MutFeature & ft) {
	if (classes.count(&ft) > 0) {
		CError error(CErrorType::InvalidArguments, "MutClassGroup::new_class", "Duplicated feature");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		MutClass * _class = new MutClass(*this, ft);
		classes[&ft] = _class; return _class;
	}
}
void MutClassGroup::add_mutant(MutClass * _class, Mutant::ID mid) {
	if (index.count(mid) > 0) {
		CError error(CErrorType::InvalidArguments, "MutClassGroup::add_mutant",
			"Duplicated mutant (" + std::to_string(mid) + ")");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		_class->add(mid); index[mid] = _class;
	}
}
void MutClassGroup::clear() {
	auto beg = classes.begin();
	auto end = classes.end();
	while (beg != end) 
		delete ((beg++)->second);

	classes.clear(); index.clear();
}
MutClass & MutClassGroup::get_class(MutFeature & ft) const {
	if (classes.count(&ft) == 0) {
		CError error(CErrorType::InvalidArguments, "MutClassGroup::get_class", "Undefined feature");
		CErrorConsumer::consume(error); exit(CErrorType::InvalidArguments);
	}
	else {
		auto iter = classes.find(&ft); return *(iter->second);
	}
}

void MutClassifierByOperator::classify_mutants() {
	/* getters */
	MutantSpace & mspace = group.get_mutant_space();
	Mutant::ID mid = 0, mnum = mspace.number_of_mutants();

	/* iterate each mutant in space */
	while (mid < mnum) {
		/* get the next mutant */
		Mutant & mutant = mspace.get_mutant(mid++);
		const std::string & oprt = mutant.get_operator();

		/* get the mutant class */
		MutClass * _class; 
		if (class_map.count(oprt) == 0) {
			MutFeature * feature = new MutFeature(oprt);
			_class = this->new_class(*feature); 
			class_map[oprt] = _class;
		}
		else {
			auto iter = class_map.find(oprt);
			_class = iter->second;
		}

		/* insert the mutant into the class */
		this->add_mutant(_class, mutant.get_id());
	} /* end while: mid < mnum */
}
void MutClassifierByLocation::classify_mutants() {
	/* getters */
	MutantSpace & mspace = group.get_mutant_space();
	Mutant::ID mid = 0, mnum = mspace.number_of_mutants();

	/* iterate each mutant in space */
	while (mid < mnum) {
		/* get the next mutant and its location (key) */
		Mutant & mutant = mspace.get_mutant(mid++);
		const Mutation & mutation = mutant.get_mutation(mutant.get_orders() - 1);
		const CodeLocation & loc = mutation.get_location();
		std::string loc_key = std::to_string(loc.get_bias());
		loc_key += ":" + std::to_string(loc.get_length());

		/* get the mutant class */
		MutClass * _class;
		if (class_map.count(loc_key) == 0) {
			MutFeature * feature = new MutFeature(loc);
			_class = this->new_class(*feature);
			class_map[loc_key] = _class;
		}
		else {
			auto iter = class_map.find(loc_key);
			_class = iter->second;
		}

		/* insert the mutant into the class */
		this->add_mutant(_class, mutant.get_id());
	} /* end while: mid < mnum */
}
void MutClassifierByCoverage::classify_mutants() {
	CoverageVector * covvec;
	while ((covvec = producer.produce()) != nullptr) {
		/* get the leaf for the coverage */
		const BitSeq & coverage = covvec->get_coverage();
		BitTrie * leaf = trie.insert_vector(coverage);

		/* get its mutant class */
		MutClass * _class;
		if (leaf->get_data() == nullptr) {
			MutFeature * feature = new MutFeature(coverage);
			_class = this->new_class(*feature);
			leaf->set_data(_class);
		}
		else _class = (MutClass *)(leaf->get_data());

		/* insert mutant */ 
		this->add_mutant(_class, covvec->get_mutant());

		consumer.consume(covvec);	/* consume the vector */
	}
}
void MutClassifierByScore::classify_mutants() {
	ScoreVector * scrvec;
	while ((scrvec = producer.produce()) != nullptr) {
		/* get the leaf for the coverage */
		const BitSeq & killset = scrvec->get_vector();
		BitTrie * leaf = trie.insert_vector(killset);

		/* get its mutant class */
		MutClass * _class;
		if (leaf->get_data() == nullptr) {
			MutFeature * feature = new MutFeature(killset);
			_class = this->new_class(*feature);
			leaf->set_data(_class);
		}
		else _class = (MutClass *)(leaf->get_data());

		/* insert mutant */
		this->add_mutant(_class, scrvec->get_mutant());

		consumer.consume(scrvec);	/* consume the vector */
	}
}

int main() { 
	// initialization
	std::string prefix = "../../../MyData/SiemensSuite/"; std::string prname = "tcas";
	File & root = *(new File(prefix + prname)); TestType ttype = TestType::tcas;

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

	// classifier 
	cfile_beg = cfiles.begin(), cfile_end = cfiles.end();
	while (cfile_beg != cfile_end) {
		// get next file and the mutant space 
		CodeFile & cfile = *(*(cfile_beg++)); cspace.load(cfile);
		MutantSpace & mspace = cmutant.get_mutants_of(cfile);

		// classifier
		MutClassGroup group(mspace);
		MutClassifierByOperator classifier(group);
		classifier.classify();

		std::cout << "Classify for \"" << cfile.get_file().get_path() << "\"\n";

		const std::map<MutFeature *, MutClass *> & classes = group.get_classes();
		auto beg = classes.begin(), end = classes.end(); size_t num = 0;
		while (beg != end) {
			MutClass & ci = *((beg++)->second); num += ci.get_mutants().number_of_mutants();
			std::cout << "\t" << ci.get_feature().get_string() << " : \t" << ci.get_mutants().number_of_mutants() << std::endl;
		}
		std::cout << "\n\tSummary: \t" << num << std::endl;
	}

	/* end */
	std::cout << "\nPress any key to exit...";
	getchar(); return 0;
}


