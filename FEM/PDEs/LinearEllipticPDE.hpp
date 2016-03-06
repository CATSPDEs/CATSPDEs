#pragma once
#include "Point.hpp"
#include "SymmetricMatrix.hpp"

class LinearEllipticPDE {
protected:
	SymmetricMatrix _diffusionTensor;
	REAL (*_forceTerm)(Point const &);
public:
	LinearEllipticPDE(SymmetricMatrix const & diffusionTensor, REAL (*forceTerm)(Point const &)) 
		: _diffusionTensor(diffusionTensor)
		, _forceTerm(forceTerm) {
	}
	SymmetricMatrix getDiffusionTensor() { return _diffusionTensor; }
	REAL getDiffusionTensor(size_t i, size_t j) { return _diffusionTensor.get(i, j); }
	REAL getForceTerm(Point const & point) { return _forceTerm(point); }
};