#pragma once
#include <iostream>
#include "Point.hpp"

typedef double(*Function)(Point const &); // function pointers

// What does it mean to read Function from stdin? 
// No sense in reading pointers, so we define our function as constant zero 

inline double zeroFunc(Point const &) { return 0.; }

inline std::istream& operator>>(std::istream& input, Function f) {
	f = zeroFunc;
	return input;
}