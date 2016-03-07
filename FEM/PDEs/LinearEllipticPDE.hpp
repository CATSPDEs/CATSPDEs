#pragma once
#include "Function.hpp"
#include "SymmetricContainer.hpp"

typedef SymmetricContainer<double> ConstTensor;
typedef SymmetricContainer<Function> Tensor;

class LinearEllipticPDE { // div(D grad u) + gamma u = f, D is constant a symmetric diffusion matrix, gamma is a constant, f is a function
	ConstTensor _diffusionTensor;
	Function _forceTerm;
	double _gamma;
public:
	LinearEllipticPDE(ConstTensor const & diffusionTensor, double gamma, Function forceTerm)
		: _diffusionTensor(diffusionTensor)
		, _gamma(gamma)
		, _forceTerm(forceTerm) {
	}
	// set / get methods
	ConstTensor diffusionTensor() const { return _diffusionTensor; }
	double diffusionTensor(size_t i, size_t j) const { return _diffusionTensor(i, j); }
	double gamma() const { return _gamma; }
	double forceTerm(Point const & point) const { return _forceTerm(point); }
};

// this is how one can implement something more general:
class NonLinearEllipticPDE {
	Tensor _diffusionTensor;
	Function _forceTerm;
	Function _gamma;
public:
	NonLinearEllipticPDE(Tensor const & diffusionTensor, Function gamma, Function forceTerm)
		: _diffusionTensor(diffusionTensor)
		, _gamma(gamma)
		, _forceTerm(forceTerm) {
	}
	ConstTensor diffusionTensor(Point const & point) {
		size_t n = _diffusionTensor.getOrder();
		ConstTensor D(n);
		for (size_t i = 0; i < n; ++i)
			for (size_t j = i; j < n; ++j)
				D(i, j) = _diffusionTensor(i, j)(point);
		return D;
	}
	double gamma(Point const & point) const { return _gamma(point); }
	double forceTerm(Point const & point) const { return _forceTerm(point); }
};