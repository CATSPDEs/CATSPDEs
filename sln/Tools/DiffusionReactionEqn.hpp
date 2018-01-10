#pragma once
#include "Mapping.hpp"

template <LocalIndex D>
class DiffusionReactionEqn { // –∇.(a ∇u) + cu = f
	ScalarField<D> _a, _c, _f;
	VectorField<D> _w;
public:
	explicit DiffusionReactionEqn(
		ScalarField<D> const & a,
		ScalarField<D> const & c,
		ScalarField<D> const & f,
		VectorField<D> const & w
	)
		: _a(a)
		, _c(c)
		, _f(f)
		, _w(w) {}
	auto& diffusionTerm() { return _a; }
	auto& reactionTerm () { return _c; }
	auto& forceTerm    () { return _f; }
	auto& convectionTerm() { return _w; }
	auto diffusionTerm (Node<D> const & p) const { return _a(p); }
	auto reactionTerm  (Node<D> const & p) const { return _c(p); }
	auto forceTerm     (Node<D> const & p) const { return _f(p); }
	auto convectionTerm(Node<D> const & p) const { return _w(p); }
};

	using DiffusionReactionEqn2D = DiffusionReactionEqn<2>;
	using DiffusionReactionEqn3D = DiffusionReactionEqn<3>;