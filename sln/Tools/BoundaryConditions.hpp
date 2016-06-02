#pragma once
#include "Function.hpp"

class BoundaryConditions { // boundary conditions 
	Function _D, _N, _R;
	// –n . (a nabla u) = _R (u – _D) – _N 
public:
	explicit BoundaryConditions(Function D = zeroFunc, Function N = zeroFunc, Function R = zeroFunc)
	// we have hom. Neumann conditions by default: 
	// n . (a nabla u) = 0 everywhere on the bndry
		: _D(D)
		, _N(N)
		, _R(R) {}

	Function DirichletCondition() const { return _D; }
	Function NeumannValue      () const { return _N; }
	Function RobinCoefficient  () const { return _R; }
	double DirichletCondition(Node& p) const { return _D(p); }
	double NeumannValue      (Node& p) const { return _N(p); }
	double RobinCoefficient  (Node& p) const { return _R(p); }
};
