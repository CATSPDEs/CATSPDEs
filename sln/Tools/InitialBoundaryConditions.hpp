#pragma once
#include "Function.hpp"
#include "TimeFunction.hpp"

class InitialBoundaryConditions { // boundary conditions 
	TimeFunction _D, _N, _R;
	Function _v; // initial velocity
	// –n[r] . (a[r] nabla u[r,t]) = _R[r,t] (u[r,t] – _D[r,t]) – _N[r,t],
	// u[r,0]   = _D[r,0],
	// u_t[r,0] = _v[r]
public:
	explicit InitialBoundaryConditions(TimeFunction D = zeroTimeFunc, TimeFunction N = zeroTimeFunc, TimeFunction R = zeroTimeFunc,
		                               Function v = zeroFunc)
		// we have hom. Neumann conditions by default: 
		// n . (a nabla u) = 0 everywhere on the bndry
		: _D(D)
		, _N(N)
		, _R(R) 
		, _v(v) {}
	TimeFunction DirichletCondition() const { return _D; }
	TimeFunction NeumannValue() const { return _N; }
	TimeFunction RobinCoefficient() const { return _R; }
	Function initialVelocity() const { return _v; }
	double initialVelocity(Node& p) const { return _v(p); }
	double DirichletCondition(Node& p, double t) const { return _D(p, t); }
	double NeumannValue(Node& p, double t) const { return _N(p, t); }
	double RobinCoefficient(Node& p, double t) const { return _R(p, t); }
};
