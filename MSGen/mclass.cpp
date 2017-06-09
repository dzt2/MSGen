#include "mclass.h"

MutClass::MutClass(MutClassGroup & grp, const std::string & desc) 
	: group(grp), description(desc) {
	mutants = grp.get_mutant_space().create_set();
}
MutClass::~MutClass() {
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
MutClass * MutClassGroup::new_class(const std::string & desc) {
	MutClass * _class = new MutClass(*this, desc);
	classes.push_back(_class); return _class;
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
	while (beg != end) delete *(beg++);

	classes.clear(); index.clear();
}





