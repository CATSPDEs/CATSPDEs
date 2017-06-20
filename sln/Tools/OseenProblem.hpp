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
	auto& massCoef() { return _alpha; }
	auto& windField() { return _w; }
	auto& inverseReynoldsNumber() { return _ReInv; }
	auto& forceTerm() { return _f; }
	auto& continuityTerm() { return _g; }
	auto massCoef() const { return _alpha; }
	auto windField(Node<D> const & p) const { return _w(p); }
	auto inverseReynoldsNumber() const { return _ReInv; }
	auto forceTerm(Node<D> const & p) const { return _f(p); }
	auto continuityTerm(Node<D> const & p) const { return _g(p); }
};

	using OseenProblem2D = OseenProblem<2>;
	using OseenProblem3D = OseenProblem<3>;