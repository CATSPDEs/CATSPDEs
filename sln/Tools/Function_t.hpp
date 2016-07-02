#pragma once
#include "Node.hpp"

typedef double(*Function_t)(Node&, double); // function pointer

inline double emptyFunc_t(Node&, double) { return 0.; } // for default parameters
inline double zeroFunc_t (Node&, double) { return 0.; }
inline double oneFunc_t  (Node&, double) { return 1.; }