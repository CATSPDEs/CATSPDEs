#pragma once
#include "AbstractPDE.hpp"
#include "SymmetricContainer.hpp"

typedef SymmetricContainer<double> ConstTensor;
typedef SymmetricContainer<Function> Tensor;

class LinearEllipticPDE : public AbstractPDE { // div(D grad u) + gamma u = f, D is constant a symmetric diffusion matrix, gamma is a constant, f is a function
	ConstTensor* _diffusionTensor;
	Function _forceTerm;
	double _gamma;
public:
	LinearEllipticPDE(ConstTensor* diffusionTensor, double gamma, Function forceTerm)
		: AbstractPDE(diffusionTensor->getOrder())
		, _diffusionTensor(diffusionTensor)
		, _gamma(gamma)
		, _forceTerm(forceTerm) {
	}
	// set / get methods
	ConstTensor diffusionTensor() const { return *_diffusionTensor; }
	double diffusionTensor(size_t i, size_t j) const { return (*_diffusionTensor)(i, j); }
	double gamma() const { return _gamma; }
	double forceTerm(Point const & point) const { return _forceTerm(point, 0.); }
};

// this is how one can implement something more general:
class NonLinearEllipticPDE {
	Tensor* _diffusionTensor;
	Function _forceTerm;
	Function _gamma;
public:
	NonLinearEllipticPDE(Tensor* diffusionTensor, Function gamma, Function forceTerm)
		: _diffusionTensor(diffusionTensor)
		, _gamma(gamma)
		, _forceTerm(forceTerm) {
	}
	ConstTensor diffusionTensor(Point const & point) {
		size_t n = (*_diffusionTensor).getOrder();
		ConstTensor D(n);
		for (size_t i = 0; i < n; ++i)
			for (size_t j = i; j < n; ++j)
				D(i, j) = (*_diffusionTensor)(i, j)(point, 0.);
		return D;
	}
	double gamma(Point const & point) const { return _gamma(point, 0.); }
	double forceTerm(Point const & point) const { return _forceTerm(point, 0.); }
};