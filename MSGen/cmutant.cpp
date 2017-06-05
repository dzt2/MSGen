#include "cmutant.h"

const Mutation & Mutant::get_mutation(size_t k) const {
	if (k >= orders) {
		CError error(CErrorType::OutOfIndex, "Mutant::get_mutation", "Invalid index: k = " + std::to_string(k));
		CErrorConsumer::consume(error);
	}
	return *(mutations[k]);
}
Mutant::Mutant(const MutantSpace & spac, Mutant::ID id, const std::string & op)
	: space(spac), mid(id), oprt(op), orders(0), mutations(nullptr) {}
void Mutant::set_orders(size_t ords) {
	clear_mutations();
	orders = ords;
	if (ords > 0) {
		mutations = new Mutation *[ords];
		for (size_t k = 0; k < ords; k++)
			mutations[k] = nullptr;
	}
}
void Mutant::set_mutation(size_t k, size_t bias, size_t length, const std::string & replace) {
	if (k >= orders) {
		CError error(CErrorType::OutOfIndex, "Mutant::get_mutation", "Invalid index: k = " + std::to_string(k));
		CErrorConsumer::consume(error);
	}
	if (mutations[k] != nullptr) delete mutations[k];
	mutations[k] = new Mutation(*this, space.get_code_file(), bias, length, replace);
}
void Mutant::clear_mutations() {
	for (size_t k = 0; k < orders; k++)
		if (mutations[k] != nullptr)
			delete mutations[k];
	delete mutations;
	mutations = nullptr;
}

MutantSource::MutantSource(const File & d) : dir(d) {
	if (!dir.is_directory()) {
		CError error(CErrorType::InvalidArguments,
			"MutantSource::MutantSource", "Not directory: " + d.get_path());
		CErrorConsumer::consume(error);
	}
	else {
		/* initialization */
		code = nullptr; schema = nullptr;

		/* derive files to code and schema */
		auto files = dir.list_files();
		auto beg = files.begin(), end = files.end();
		while (beg != end) {
			/* get next file */
			const File & file = *(*(beg++));
			if (file.is_directory()) continue;

			/* check its name and set code|schema */
			const std::string & name = file.get_local_name();
			if (endswith(name, ".c") && code == nullptr)
				code = &file;
			else if (name == "mschema.txt" && schema == nullptr)
				schema = &file;
		}

		/* validation */
		std::string error;
		if (code == nullptr) error = "__" + dir.get_local_name() + ".c";
		else if (schema == nullptr) error = "mschema.txt";
		else return;

		CError errors(CErrorType::Runtime, "MutantSource::MutantSource", error + " is not found");
		CErrorConsumer::consume(errors);
	}
}

MutantSpace::MutantSpace(const CMutant & p, const File & dir, const CodeFile & cfile)
	: project(p), source(dir), code_file(cfile), mutants(), mutsets() {}
Mutant * MutantSpace::create_mutant(const std::string & oprt) {
	Mutant * mutant = new Mutant(*this, mutants.size(), oprt);
	mutants.push_back(mutant); return mutant;
}
void MutantSpace::clear_mutants() {
	auto beg = mutants.begin();
	auto end = mutants.end();
	while (beg != end) {
		delete *(beg++);
	}
	mutants.clear();
}
void MutantSpace::clear_sets() {
	auto beg = mutsets.begin();
	auto end = mutsets.end();
	while (beg != end) {
		delete *(beg++);
	}
	mutsets.clear();
}
Mutant & MutantSpace::get_mutant(Mutant::ID mid) const {
	if (mid >= mutants.size()) {
		CError error(CErrorType::InvalidArguments, "MutantSpace::get_mutant", "Out of index: mid = " + std::to_string(mid));
		CErrorConsumer::consume(error);
	}
	return *(mutants[mid]);
}
MutantSet * MutantSpace::create_set() {
	MutantSet * mset = new MutantSet(*this);
	mutsets.insert(mset); return mset;
}
bool MutantSpace::delete_set(MutantSet * mset) {
	if (mutsets.count(mset) == 0) return false;
	else {
		delete mset; mutsets.erase(mset); return true;
	}
}

bool MutantSet::has_mutant(Mutant::ID mid) const {
	if (mid >= vec.bit_number())
		return false;
	else return vec.get_bit(mid) == BIT_1;
}
bool MutantSet::add_mutant(Mutant::ID mid) {
	if (vec.get_bit(mid) == BIT_1) return false;
	else {
		vec.set_bit(mid, BIT_1); number++;
		return true;
	}
}
bool MutantSet::del_mutant(Mutant::ID mid) {
	if (vec.get_bit(mid) == BIT_0)
		return false;
	else {
		number--; vec.set_bit(mid, BIT_0);
		return true;
	}
}
bool MutantSet::complement() {
	byte * bytes = vec.get_bytes();
	int bnum = vec.byte_number();

	for (int i = 0; i < bnum; i++) {
		bytes[i] = ~bytes[i];
	}

	number = vec.bit_number() - number;
	return true;
}

CMutant::CMutant(const File & root, const CodeSpace & cs) : cspace(cs), loader() {
	/* validation */
	if (!root.is_directory()) {
		CError error(CErrorType::InvalidArguments, "CMutant::CMutant", "Not directory: " + root.get_path());
		CErrorConsumer::consume(error);
	}

	/* get ../muta/ */
	dir = nullptr; auto children = root.list_files();
	auto beg = children.begin(), end = children.end();
	while (beg != end) {
		const File * file = *(beg++);
		if (!(file->is_directory())) continue;
		const std::string & name = file->get_local_name();

		if (name == "muta") {
			dir = file; break;
		}
	}
	if (dir == nullptr) {
		CError error(CErrorType::InvalidArguments, "CMutant::CMutant", "/muta is not found!");
		CErrorConsumer::consume(error);
	}

	/* find valid code-files */
	std::map<std::string, const CodeFile *> name_code;
	auto cfiles = cspace.get_code_set();
	auto cbeg = cfiles.begin(), cend = cfiles.end();
	while (cbeg != cend) {
		const CodeFile * cfile = *(cbeg++);
		CodeFileType type = cfile->get_type();
		if (type == CodeFileType::Source || type
			== CodeFileType::Preprocessed) {
			std::string name = cfile->get_file().get_local_name();
			name = name.substr(0, name.length() - 2);
			name_code[name] = cfile;
		}
	}

	/* create mutant space */
	children = dir->list_files();
	beg = children.begin(), end = children.end();
	while (beg != end) {
		const File * ndir = *(beg++);
		if (!(ndir->is_directory())) continue;
		const std::string & dname = ndir->get_local_name();

		if (name_code.count(dname) > 0) {
			auto iter = name_code.find(dname);
			const CodeFile * cfile = iter->second;

			MutantSpace * space = new MutantSpace(*this, *ndir, *cfile);
			spaces[cfile] = space;
		}
	}
	name_code.clear();
}
CMutant::~CMutant() {
	auto beg = spaces.begin(), end = spaces.end();
	while (beg != end) {
		MutantSpace * space = (beg++)->second;
		delete space;
	}
	spaces.clear();
}

void MutantLoader::load(MutantSpace & space, bool deep) const {
	/* declarations */
	const File & schema = space.get_source().get_schema();
	LineReader reader(schema.get_path()); Mutant::ID mid;
	size_t k, orders, bias, length; std::string replace, oprt;

	/* read mutant and put into space */
	while (reader.hasNext()) {
		if (get_id(reader, mid)) {	/* get next mutant id */
			if (get_order(reader, orders)) {	/* get the mutant orders */
				if (get_operator(reader, oprt)) {	/* get mutant operator */
					Mutant * mutant = space.create_mutant(oprt);
					mutant->set_orders(orders);
					if (deep) {	/* deep copy: load mutations */
						k = 0;
						while (k < orders) {
							if (get_mutations(reader, bias, length, replace)) {	/* get next mutations for mutant */
								mutant->set_mutation(k, bias - 1, length, replace);
							}
							else {
								const Mutation & mutation = mutant->get_mutation(k - 1);
								bias = mutation.get_location().get_bias();
								length = mutation.get_location().get_length();
								replace = mutation.get_replacement();
								mutant->set_mutation(k, bias, length, replace);
							}
							k = k + 1;
						}
					}	/* end if: deep copy */
				}
			}
		}
	} /* end while */
}
bool MutantLoader::get_id(LineReader & reader, Mutant::ID & mid) const {
	/* initialization */
	mid = -1; std::string linestr;

	/* read each line to locate next "MUTANT # ID" */
	while (reader.hasNext()) {
		/* get next line */
		linestr = reader.next();
		trimstring(linestr);

		if (linestr.empty()) continue;
		else if (startswith(linestr, "MUTANT #")) {	/* find next id */
			linestr = linestr.substr(8, linestr.length() - 8);
			mid = std::stoul(linestr); return true;
		}
	}

	/* not found */ return false;
}
bool MutantLoader::get_order(LineReader & reader, size_t & orders) const {
	/* initialization */
	orders = 0; std::string linestr;
	int s, num;

	/* read each line to locate next "Calling(Called) function starts at:"*/
	while (reader.hasNext()) {
		linestr = reader.next();	/* get next line */
		trimstring(linestr);

		/* locate target line */
		if (startswith(linestr, "Calling function starts at:")
			|| startswith(linestr, "Called function starts at:")) {
			/* get function bias */
			s = index_of(linestr, ':') + 1;
			linestr = linestr.substr(s, linestr.length() - s);
			num = std::stoi(linestr);

			if (num >= 0) orders++;	/* -1 means no fault is seeded at the location */
		}
		else if (startswith(linestr, "Last test case used")) break;
	}

	/* found */ return orders > 0;
}
bool MutantLoader::get_operator(LineReader & reader, std::string & oprt) const {
	oprt = ""; std::string linestr; int s, e;

	while (reader.hasNext()) {
		linestr = reader.next();
		trimstring(linestr);

		if (startswith(linestr, "Operator:")) {
			s = index_of(linestr, '(');
			e = index_of(linestr, ')');
			oprt = linestr.substr(s + 1, e - s - 1);
			return true;
		}
		else if (startswith(linestr, "Descriptor:")
			|| linestr.empty()) {
			break;
		}
	}

	return false;
}
bool MutantLoader::get_mutations(LineReader & reader,
	size_t & bias, size_t & length, std::string & replace) const {
	/* initialization */
	bias = 0, length = 0, replace = "";
	std::string linestr, numstr; int s, e;

	/* iterate each line to get next mutation */
	while (reader.hasNext()) {
		linestr = reader.next();
		trimstring(linestr);

		if (startswith(linestr, "Offset:")) {
			/* get offset */
			s = index_of(linestr, ':');
			e = index_of(linestr, ',');
			numstr = linestr.substr(s + 1, e - s - 1);
			bias = std::stoul(numstr);

			/* get length */
			linestr = linestr.substr(e + 1, linestr.length() - e - 1);
			linestr = linestr.substr(9, linestr.length() - 9);
			e = index_of(linestr, 'c');
			numstr = linestr.substr(0, e);
			length = std::stoul(numstr);
		}
		else if (startswith(linestr, "Get on:")) {
			/* get replace */
			s = index_of(linestr, ':') + 1;
			replace = linestr.substr(s, linestr.length() - s);
			return true;
		}
		else if (linestr.empty()) break;
	}

	/* return */ return false;
}

/*
int main() {
// initialization
File * dir = new File("print_tokens");
CProgram * codprog = new CProgram(*dir);
CMutant * mutprog = new CMutant(*dir, codprog->get_source());

// load mutants
auto spaces = mutprog->get_spaces();
auto beg = spaces.begin(), end = spaces.end();
while (beg != end) {
// get next space
const CodeFile & cfile = *(beg->first);
MutantSpace & space = *((beg++)->second);

// before load
std::cout << cfile.get_file().get_path() <<
" : " << space.number_of_mutants() << std::endl;

// load mutants
mutprog->load_mutants_for(space, true);

// after load
std::cout << cfile.get_file().get_path() <<
" : " << space.number_of_mutants() << std::endl;

// to the next
std::cout << std::endl;
}


// release resources
delete mutprog; delete codprog; delete dir;
getchar(); return 0;
}
*/


