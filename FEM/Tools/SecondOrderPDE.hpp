#pragma once
#include "Function.hpp"

class SecondOpderPDE { // div ( lambda grad u ) + gamma u + f = 0
	Function _lambda, _gamma, _f;
public:
	explicit SecondOpderPDE(Function lambda = oneFunc, Function gamma = zeroFunc, Function f = zeroFunc) // Laplace's eqn by default
		: _lambda(lambda)
		, _gamma(gamma)
		, _f(f) {}
	virtual ~SecondOpderPDE() {}
	double lambda(Point const & p) const { return _lambda(p); }
	double gamma(Point const & p) const { return _gamma(p); }
	double f(Point const & p) const { return _f(p); }
};