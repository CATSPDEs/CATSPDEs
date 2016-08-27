#pragma once
#include <iostream>
#include <ctime>
#include <string>
#include <stack>
#include "rlutil.h" // https://github.com/tapio/rlutil
// for cross–platform terminal text coloring

class SingletonLogger {
	std::stack<time_t> _processes;
	double _diff(clock_t t1, clock_t t2) const { return fabs(t2 - t1) / CLOCKS_PER_SEC; }
	// singleton
	SingletonLogger();
	SingletonLogger(SingletonLogger const &);
	SingletonLogger& operator=(SingletonLogger const &);
public:
	~SingletonLogger();
	static SingletonLogger& instance();
	void beg(std::string const &); // begin new process
	void end(); // terminate last process
	void log(std::string const &) const; // log message
	void mes(std::string const &) const; // log message (just tabs)
	void err(std::string const &) const; // print error
	void inp(std::string const &) const; // print input invite
	std::string tab() const { // tabulations
		return std::string(_processes.size(), '\t'); 
	} 
};