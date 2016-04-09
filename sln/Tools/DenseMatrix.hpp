#pragma once
#include "AbstractSquareMatrix.hpp"
#include "IRealMatrix.hpp"

class DenseMatrix : public AbstractSquareMatrix<double>, public IRealMatrix {
	double** _A, * _beg;
	// to be implemented
	double& _set(size_t, size_t);
	double _get(size_t, size_t) const;
public:
	explicit DenseMatrix(size_t n = 0);
	~DenseMatrix();
	// to be implemented
	std::vector<double> solve(std::vector<double> const &);
	std::vector<double> mult(std::vector<double> const &) const;
};
