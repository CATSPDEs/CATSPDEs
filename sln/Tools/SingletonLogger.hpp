#pragma once
#include <iostream>
#include <sstream>
#include <ctime>
#include <string>
#include <stack>
#include <vector>
#include "rlutil.h" // https://github.com/tapio/rlutil
// for cross–platform terminal text coloring

class SingletonLogger {
	std::stack<time_t> _processes; // stack of started proccesses
	std::stringstream _stream;
	double _diff(clock_t t1, clock_t t2) const { return fabs(t2 - t1) / CLOCKS_PER_SEC; }
	std::string _format(std::string const &) const;
	// singleton
	SingletonLogger();
	SingletonLogger(SingletonLogger const &);
	SingletonLogger& operator=(SingletonLogger const &);
public:
	std::ostringstream buf;
	~SingletonLogger();
	static SingletonLogger& instance();
	void beg(std::string const &); // begin new process
	void end(); // terminate last process
	void log(std::string const &) const; // print message w/ "[log]" prefix
	void log(); // print formatted buf 
	void wrn(std::string const &) const; // print error
	void err(std::string const &) const; // print error
	void inp(std::string const &) const; // print input invite
	bool yes(std::string const &) const; // true if get 'y' from stdin
	size_t opt(std::string const &, std::vector<std::string> const &); // choose vector element (get index from stdin), return its index
	std::string tab() const { // tabulations
		return std::string(_processes.size(), '\t'); 
	} 
};