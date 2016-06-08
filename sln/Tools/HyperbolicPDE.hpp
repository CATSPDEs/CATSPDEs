#pragma once
#include "Function.hpp"
#include "TimeFunction.hpp"

class HyperbolicPDE { 
	// chi[r] u_tt[r,t] + sigma[r] u_t[r,t] – nabla . (a[r] nabla u[r,t]) = f[r,t],
	// where r := <x, y> — space variable, t := time, and subscripts mean derivatives
	Function _chi, _sigma, _a;
	TimeFunction _f;
public:
	explicit HyperbolicPDE(Function chi = oneFunc,
	                       Function sigma = zeroFunc,
	                       Function a = oneFunc, 
	                       TimeFunction f = zeroTimeFunc)
		// hom. wave eqn by default
		: _chi(chi)
		, _sigma(sigma)
		, _a(a)
		, _f(f) {}
	virtual ~HyperbolicPDE() {}
	Function chi() const { return _chi; }
	Function sigma() const { return _sigma; }
	Function diffusionTerm() const { return _a; }
	TimeFunction forceTerm() const { return _f; }
	double chi(Node& p) const { return _chi(p); }
	double sigma(Node& p) const { return _sigma(p); }
	double diffusionTerm(Node& p) const { return _a(p); }
	double forceTerm(Node& p, double t) const { return _f(p, t); }
};