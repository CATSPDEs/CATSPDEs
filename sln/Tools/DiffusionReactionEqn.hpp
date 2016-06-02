#pragma once
#include "Function.hpp"

class DiffusionReactionEqn { // –nabla . (a nabla u) + cu = f
	Function _a, _c, _f;
public:
	explicit DiffusionReactionEqn(Function a = oneFunc, Function c = zeroFunc, Function f = zeroFunc) 
		// Laplace’s eqn by default
		: _a(a)
		, _c(c)
		, _f(f) {}
	virtual ~DiffusionReactionEqn() {}
	Function diffusionTerm() const { return _a; }
	Function reactionTerm () const { return _c; }
	Function forceTerm    () const { return _f; }
	double diffusionTerm(Node& p) const { return _a(p); }
	double reactionTerm (Node& p) const { return _c(p); }
	double forceTerm    (Node& p) const { return _f(p); }
};