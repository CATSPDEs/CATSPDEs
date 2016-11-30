#pragma once
#include "Mapping.hpp"

template <LocalIndex D>
class OseenProblem { // alpha u + (w.∇)u – Re^-1 Δu + ∇p = f, ∇.u = g
	double _alpha, _ReInv;
	VectorField<D> _w, _f;
	ScalarField<D> _g;
public:
	explicit OseenProblem( // hom. Stokes problem by default
		double alpha,
		VectorField<D> const & w,
		double ReInv,
		VectorField<D> const & f,
		ScalarField<D> const & g
	) 
		: _alpha(alpha)
		, _w(w)
		, _ReInv(ReInv)
		, _f(f)
		, _g(g) {}
	auto massTerm() const { return _alpha; }
	auto NewtonTerm() const { return _w; }
	auto inverseReynoldsNumber() const { return _ReInv; }
	auto ReynoldsNumber() const { return 1. / _ReInv; }
	auto forceTerm() const { return _f; }
	auto continuityTerm() const { return _g; }
	auto NewtonTerm(Node<D> const & p) const { return _w(p); }
	auto forceTerm(Node<D> const & p) const { return _f(p); }
	auto continuityTerm(Node<D> const & p) const { return _g(p); }
};

	using OseenProblem2D = OseenProblem<2>;
	using OseenProblem3D = OseenProblem<3>;