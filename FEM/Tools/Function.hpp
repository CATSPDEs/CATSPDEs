#pragma once
#include <iostream>
#include <stdexcept>
#include "Point.hpp"

typedef double(*Function)(Point const &); // function pointer

inline double zeroFunc(Point const &) { return 0.; }
inline double oneFunc(Point const &) { return 1.; }

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