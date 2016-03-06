#pragma once
#include "Point.hpp"
#include "SymmetricMatrixOfFunctions.hpp"

class NonLinearEllipticPDE {
	SymmetricMatrixOfFunctions _diffusionTensor;
	fnctn _forceTerm;
	fnctn _gamma;
public:
	NonLinearEllipticPDE(SymmetricMatrixOfFunctions const & diffusionTensor, fnctn gamma, fnctn forceTerm)
		: _diffusionTensor(diffusionTensor)
		, _gamma(gamma)
		, _forceTerm(forceTerm) {
	}
	double getDiffusionTensor(size_t i, size_t j, Point const & point) { return _diffusionTensor(i, j)(point); }
	double getGamma(Point const & point) { return _gamma(point); }
	double getForceTerm(Point const & point) { return _forceTerm(point); }
};