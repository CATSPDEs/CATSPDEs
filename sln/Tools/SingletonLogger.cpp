#include "SingletonLogger.hpp"

SingletonLogger::SingletonLogger() {
	rlutil::saveDefaultColor();
	rlutil::setBackgroundColor(0);
}

SingletonLogger::~SingletonLogger() {
	rlutil::resetColor();
}

SingletonLogger& SingletonLogger::instance() {
	static SingletonLogger single;
	return single;
}

void SingletonLogger::beg(std::string const & message) {
	rlutil::setColor(11);
	std::cout  << tab() << "[beg] ";
	rlutil::setColor(15);
	std::cout << message << " . . .\n";
	rlutil::setColor(8);
	_processes.push(clock());
}

void SingletonLogger::end() {
	if (_processes.size() == 0) {
		err("there is nothing to end");
		return;
	}
	clock_t t = _processes.top();
	_processes.pop();
	rlutil::setColor(11);
	std::cout << tab() << "[end] ";
	rlutil::setColor(15);
	std::cout << _diff(clock(), t) << '\n';
	rlutil::setColor(8);
}

void SingletonLogger::err(std::string const & message) const {
	rlutil::setColor(12);
	std::cout << tab() << "[err] ";
	rlutil::setColor(8);
	std::cout << message << '\n';
}

void SingletonLogger::mes(std::string const & message) const {
	std::cout << tab() << "      " << message << '\n';
}

void SingletonLogger::log(std::string const & message) const {
	rlutil::setColor(13);
	std::cout << tab() << "[log] ";
	rlutil::setColor(8);
	std::cout << message << '\n';
}

void SingletonLogger::inp(std::string const & message) const {
	rlutil::setColor(14);
	std::cout << tab() << "[inp] ";
	rlutil::setColor(8);
	std::cout << message << ":\n" << tab();
}