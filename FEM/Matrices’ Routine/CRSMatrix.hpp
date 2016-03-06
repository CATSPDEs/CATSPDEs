// declaration of CRSMatrix

#pragma once
#include "AbstractRealMatrix.hpp"

class CRSMatrix : public AbstractRealMatrix { // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	std::vector<size_t> _ptr, // ptr := locations in the val vector that start a row
						_col; // col := column indices of the elements in the val vector
	std::vector<double> _val; // val := values of the nonzero elements of the matrix (row-wised)
public:
	CRSMatrix(size_t, size_t); // order and number of nonzeros
	~CRSMatrix() {}
	// virtual methods to be implemented
	double& set(size_t, size_t);
	double get(size_t, size_t) const;
	std::vector<double> solve(std::vector<double> const &);
	std::istream& load(std::istream&);
	std::ostream& save(std::ostream&) const;
	std::vector<double> mult(std::vector<double> const &) const;
	// new methods
	std::vector<double> CG(std::vector<double> const &, double e = 1E-14, unsigned count = 10000U) const; // conjugate gradients
	friend CRSMatrix operator*(CRSMatrix const &, CRSMatrix const &); // matrix-matrix product
};