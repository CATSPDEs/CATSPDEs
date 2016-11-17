#pragma once
#include "Mapping.hpp"

template <LocalIndex D>
class DiffusionReactionEqn { // –∇.(a ∇u) + cu = f
	ScalarField<D> _a, _c, _f;
public:
	explicit DiffusionReactionEqn(ScalarField<D> a = constUnity<D, 1>, ScalarField<D> c = constZero<D, 1>, ScalarField<D> f = constZero<D, 1>)
		// Laplace’s eqn by default
		: _a(a)
		, _c(c)
		, _f(f) {}
	~DiffusionReactionEqn() {}
	auto diffusionTerm() const { return _a; }
	auto reactionTerm () const { return _c; }
	auto forceTerm    () const { return _f; }
	auto diffusionTerm(Node<D> const & p) const { return _a(p); }
	auto reactionTerm (Node<D> const & p) const { return _c(p); }
	auto forceTerm    (Node<D> const & p) const { return _f(p); }
};

using DiffusionReactionEqn2D = DiffusionReactionEqn<2>;