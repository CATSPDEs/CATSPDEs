#pragma once
#include "Node.hpp"

typedef double(*TimeFunction)(Node&, double); // function pointer

inline double zeroTimeFunc(Node&, double) { return 0.; }
inline double oneTimeFunc(Node&, double) { return 1.; }