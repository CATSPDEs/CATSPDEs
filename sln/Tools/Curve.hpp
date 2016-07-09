#pragma once
#include "Node.hpp"
#include "constants.hpp"

typedef Node(*Curve)(double); // function pointer

// unit circle curve w/ center at origo
inline Node circleCurve(double theta) { // unit speed counter-clock parametrization
	return Node(cos(2 * PI * theta), sin(2 * PI * theta));
}