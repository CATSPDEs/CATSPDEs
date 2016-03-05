// declaration of CRSMatrix

#pragma once
#include "AbstractSquareMatrix.hpp"

class CRSMatrix : public AbstractSquareMatrix { // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	std::vector<size_t> _ptr, // ptr := locations in the val vector that start a row
						_col; // col := column indices of the elements in the val vector
	std::vector<REAL> _val; // val := values of the nonzero elements of the matrix (row-wised)
public:
	CRSMatrix(size_t, size_t); // order and number of nonzeros
	~CRSMatrix() {}
	// virtual methods to be implemented
	REAL& set(size_t, size_t);
	REAL get(size_t, size_t) const;
	std::vector<REAL> solve(std::vector<REAL> const &);
	std::istream& load(std::istream&);
	std::ostream& save(std::ostream&) const;
	std::vector<REAL> mult(std::vector<REAL> const &) const;
	// new methods
	std::vector<REAL> CG(std::vector<REAL> const &, REAL e = 1E-14, unsigned count = 10000U) const; // conjugate gradients
	friend CRSMatrix operator*(CRSMatrix const &, CRSMatrix const &); // matrix-matrix product
};