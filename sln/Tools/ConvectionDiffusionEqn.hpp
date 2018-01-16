#pragma once
#include "Mapping.hpp"

template <LocalIndex D>
class ConvectionDiffusionEqn { // –∇.(diffusion ∇u) + convection.∇u + reactiom u = force
	ScalarField<D> _diffusion, _reaction, _force;
	VectorField<D> _convection;
public:
	explicit ConvectionDiffusionEqn( // Laplace equation by default
		ScalarField<D> const & diffusion = unityMapping<D, 1>,
		ScalarField<D> const & reaction = zeroMapping<D, 1>,
		ScalarField<D> const & force = zeroMapping<D, 1>,
		VectorField<D> const & convection = zeroMapping<D, D>
	)
		: _diffusion(diffusion)
		, _reaction(reaction)
		, _force(force)
		, _convection(convection) {}
	auto& diffusion() { return _diffusion; }
	auto& reaction() { return _reaction; }
	auto& force() { return _force; }
	auto& convection() { return _convection; }
	auto diffusion(Node<D> const & p) const { return _diffusion(p); }
	auto reaction(Node<D> const & p) const { return _reaction(p); }
	auto force(Node<D> const & p) const { return _force(p); }
	auto convection(Node<D> const & p) const { return _convection(p); }
};

	using ConvectionDiffusionEqn1D = ConvectionDiffusionEqn<1>;
	using ConvectionDiffusionEqn2D = ConvectionDiffusionEqn<2>;
	using ConvectionDiffusionEqn3D = ConvectionDiffusionEqn<3>;