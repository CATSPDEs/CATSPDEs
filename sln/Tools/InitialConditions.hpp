#pragma once
#include "Function.hpp"

class InitialConditions { // boundary conditions 
	Function _u0; // initial position
	Function _v0; // initial velocity
	// u  [r, t0] = _u0[r],
	// u_t[r, t0] = _v0[r]
public:
	explicit InitialConditions(Function u0 = zeroFunc,
		                       Function v0 = zeroFunc)
		: _u0(u0)
		, _v0(v0) {
	}
	double initialPosition(Node& p) const { return _u0(p); }
	double initialVelocity(Node& p) const { return _v0(p); }
};
