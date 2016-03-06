#pragma once
#include "Point.hpp"
#include "SymmetricMatrix.hpp"

typedef double(*fnctn)(Point const &);

class LinearEllipticPDE { // div(D grad u) + gamma u = f, D is constant a symmetric diffusion matrix, gamma is a constant, f is a function
	SymmetricMatrix _diffusionTensor;
	fnctn _forceTerm;
	double _gamma;
public:
	LinearEllipticPDE(SymmetricMatrix const & diffusionTensor, double gamma, fnctn forceTerm)
		: _diffusionTensor(diffusionTensor)
		, _gamma(gamma)
		, _forceTerm(forceTerm) {
	}
	SymmetricMatrix getDiffusionTensor() { return _diffusionTensor; }
	double getDiffusionTensor(size_t i, size_t j) { return _diffusionTensor.get(i, j); }
	double getGamma() { return _gamma; }
	double getForceTerm(Point const & point) { return _forceTerm(point); }
};