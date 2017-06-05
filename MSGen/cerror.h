#pragma once

#include <iostream>

/* type of error in system */
enum CErrorType {
	NoErrors,
	InvalidArguments,
	OutOfIndex,
	Runtime,
};
/* object for error */
class CError {
public:
	/* create an error object */
	CError(enum CErrorType tp, const std::string & md, const std::string & rs)
		: type(tp), method(md), reason(rs) {}
	/* deconstructor */
	~CError() {}

	/* get error type */
	enum CErrorType get_type() const { return type; }
	/* get the method name, where error occurs */
	const std::string & get_method() const { return method; }
	/* get the reason why error is caused */
	const std::string & get_reason() const { return reason; }

protected:
	enum CErrorType type;
	std::string method;
	std::string reason;
};
/* consumer for errors */
class CErrorConsumer {
public:
	static void consume(const CError & err) { consume(err, std::cerr); }
	static void consume(const CError & err, std::ostream & out) {
		out << "Error {" << err.get_type() << "} at \"";
		out << err.get_method().c_str() << "\" for:\n\t"
			<< err.get_reason().c_str() << std::endl;
		throw err.get_type();
	}
};
