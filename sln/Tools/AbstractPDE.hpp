#pragma once
#include "Function.hpp"

class AbstractPDE {
	unsigned _dim;
public:
	explicit AbstractPDE(unsigned dim = 1U) : _dim(dim) {
		if (dim < 1 || 3 < dim) throw std::out_of_range("we are here to solve 1D / 2D / 3D problems");
	}
	virtual ~AbstractPDE() {}
	unsigned getDim() { return _dim; }
};
