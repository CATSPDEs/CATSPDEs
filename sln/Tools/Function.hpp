#pragma once
#include <iostream>
#include <stdexcept>
#include "Node.hpp"

typedef double(*Function)(Node&); // function pointer

inline double emptyFunc(Node&) { return 0.; } // for default parameters
inline double zeroFunc (Node&) { return 0.; }
inline double oneFunc  (Node&) { return 1.; }

// What does it mean to read Function from stdin? 
// No sense in saving pointers

inline std::istream& operator>>(std::istream& input, Function& f) {	
	unsigned n;
	input >> n;
	if (n == 0U) f = zeroFunc;
	else if (n == 1U) f = oneFunc;
	else throw std::logic_error("you want to read something wierd");
	return input;
}

inline std::ostream& operator<<(std::ostream& output, Function const & f) {
	if (f == zeroFunc) output << 0;
	else if (f == oneFunc) output << 1;
	else throw std::logic_error("you want to save something wierd");
	return output;
}