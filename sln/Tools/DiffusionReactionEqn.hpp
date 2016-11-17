#pragma once
#include "Function.hpp"

template <LocalIndex D>
class DiffusionReactionEqn { // –∇.(a ∇u) + cu = f
	Function<D> _a, _c, _f;
public:
	explicit DiffusionReactionEqn(Function<D> a = oneFunc, Function<D> c = zeroFunc, Function<D> f = zeroFunc) 
		// Laplace’s eqn by default
		: _a(a)
		, _c(c)
		, _f(f) {}
	~DiffusionReactionEqn() {}
	Function<D> diffusionTerm() const { return _a; }
	Function<D> reactionTerm () const { return _c; }
	Function<D> forceTerm    () const { return _f; }
	double diffusionTerm(Node<D> const & p) const { return _a(p); }
	double reactionTerm (Node<D> const & p) const { return _c(p); }
	double forceTerm    (Node<D> const & p) const { return _f(p); }
};

using DiffusionReactionEqn2D = DiffusionReactionEqn<2>;