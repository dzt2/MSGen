#include "cfile.h"
#include <queue>

ExeSource::ExeSource(const File & dir) : root(dir), program(nullptr), trace(nullptr) {
	if (!dir.exist() || !dir.is_directory()) {
		CError error(CErrorType::InvalidArguments, "ExeSource", "Undefined file");
		CErrorConsumer::consume(error);
	}
	else {
		auto files = dir.list_files();
		auto beg = files.begin(), end = files.end();
		while (beg != end) {
			File * file = *(beg++);
			const std::string & name = file->get_local_name();

			if (endswith(name, ".trace.exe")) {
				if (trace == nullptr) {
					trace = file;
				}
			}
			else if (endswith(name, ".exe")) {
				if (program == nullptr) {
					program = file;
				}
			}
		}

		if (program == nullptr) {
			CError error(CErrorType::InvalidArguments, "ExeSource", "No xxx.exe is found");
			CErrorConsumer::consume(error);
		}
		else if (trace == nullptr) {
			CError error(CErrorType::InvalidArguments, "ExeSource", "No xxx.trace.exe is found");
			CErrorConsumer::consume(error);
		}
	}
}

ExeFile::ExeFile(const ExeSpace & spac, const File & fil) : space(spac), file(fil) {
	if (!file.exist() || file.is_directory()) {
		CError error(CErrorType::InvalidArguments, "ExeFile", "Invalid file");
		CErrorConsumer::consume(error);
	}
}

ExeSpace::ExeSpace(const CProgram & p, const ExeSource & src) : project(p), source(src) {
	exec_program = new ExeFile(*this, src.get_exe());
	trace_program = new ExeFile(*this, src.get_trace_exe());
}

CodeSource::CodeSource(const File & d) : dir(d) {
	if (!dir.exist() || !dir.is_directory()) {
		CError error(CErrorType::InvalidArguments, "CodeSource", "Invalid directory");
		CErrorConsumer::consume(error);
	}
}

CodeFile::CodeFile(const CodeSpace & spac, const File & fil) : space(spac), source(fil), text(nullptr) {
	const std::string & name = source.get_local_name();
	if (endswith(name, ".h")) type = CodeFileType::Header;
	else if (endswith(name, ".c")) type = CodeFileType::Source;
	else if (endswith(name, ".i")) type = CodeFileType::Preprocessed;
	else {
		CError error(CErrorType::InvalidArguments, "CodeFile", "Not a code file");
		CErrorConsumer::consume(error);
	}
}

CodeSpace::CodeSpace(const CProgram & p, const CodeSource & src)
	: project(p), source(src), code_set(), loader() {
	std::queue<const File *> fqueue;
	fqueue.push(&(source.get_directory()));

	const File * file;
	while (!fqueue.empty()) {
		file = fqueue.front();
		fqueue.pop();

		if (file->is_directory()) {
			auto files = file->list_files();
			auto beg = files.begin(), end = files.end();
			while (beg != end) {
				fqueue.push(*(beg++));
			}
		}
		else {
			const std::string & name = file->get_local_name();
			if (endswith(name, ".h") || endswith(name, ".c") || endswith(name, ".i"))
				code_set.insert(new CodeFile(*this, *file));
		}
	}
}
CodeSpace::~CodeSpace() {
	auto beg = code_set.begin(), end = code_set.end();
	while (beg != end) {
		delete *(beg++);
	}
	code_set.clear();
}

void CodeLoader::load_in(CodeFile & file) const {
	if (file.text == nullptr) {
		LineReader reader(file.get_file().get_path());
		file.text = new TextBuild(reader);
	}
}
void CodeLoader::release(CodeFile & file) const {
	if (file.text != nullptr) {
		delete file.text;
		file.text = nullptr;
	}
}

CProgram::CProgram(const File & dir) : root(dir),
exe_source(nullptr), code_source(nullptr), orig_source(nullptr) {
	auto files = dir.list_files();
	auto beg = files.begin(), end = files.end();
	while (beg != end) {
		File * file = *(beg++);
		if (file->is_directory()) {
			const std::string & name = file->get_local_name();

			if (name == "exec" && exe_source == nullptr)
				exe_source = new ExeSource(*file);
			else if (name == "source" && code_source == nullptr)
				code_source = new CodeSource(*file);
			else if (name == "source.orig" && orig_source == nullptr)
				orig_source = new CodeSource(*file);
		}
	}

	if (exe_source == nullptr) {
		CError error(CErrorType::InvalidArguments, "CProgram", "/exec is not found");
		CErrorConsumer::consume(error);
	}
	else if (code_source == nullptr) {
		CError error(CErrorType::InvalidArguments, "CProgram", "/source is not found");
		CErrorConsumer::consume(error);
	}
	else if (orig_source == nullptr) {
		CError error(CErrorType::InvalidArguments, "CProgram", "/source.orig is not found");
		CErrorConsumer::consume(error);
	}

	exe_space = new ExeSpace(*this, *exe_source);
	source = new CodeSpace(*this, *code_source);
	source_orig = new CodeSpace(*this, *orig_source);
}
CProgram::~CProgram() {
	delete exe_space;
	delete source;
	delete source_orig;

	delete exe_source;
	delete code_source;
	delete orig_source;
}

std::string CodeLocation::get_text_at() const {
	if (file.get_text() == nullptr) {
		CError error(CErrorType::InvalidArguments, "CodeLocation", "Text is not loaded yet");
		CErrorConsumer::consume(error);
	}
	else {
		const TextBuild & text = *(file.get_text());
		return text.getText().substr(bias, length);
	}
}