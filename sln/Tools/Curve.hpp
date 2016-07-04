#pragma once
#include "Node.hpp"

typedef Node(*Curve)(double); // function pointer

double const PI = 3.141592653589793;

// unit circle curve w/ center at origo
inline Node circleCurve(double theta) { // unit speed counter-clock parametrization
	return Node(cos(2 * PI * theta), sin(2 * PI * theta));
}