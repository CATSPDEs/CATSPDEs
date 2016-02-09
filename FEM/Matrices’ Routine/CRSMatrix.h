// declaration of CRSMatrix

#pragma once
#include <iostream>
#include <cstddef> // size_t
#include "real.h"
#include "vector.h"

class CRSMatrix { // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	size_t _n; // n := order of square matrix
	std::vector<size_t> _ptr, // ptr := locations in the val vector that start a row
						_col; // col := column indices of the elements in the val vector
	std::vector<REAL> _val; // val := values of the nonzero elements of the matrix (row-wised)
public:
	CRSMatrix(size_t, size_t); // order and number of nonzeros
	~CRSMatrix();
	size_t getOrder() const;
	REAL operator()(size_t, size_t) const; // get value of element using indicies of traditional form of matrix
	REAL& operator()(size_t, size_t); // set —‘‘—
	std::vector<REAL> CG(std::vector<REAL> const &, REAL e = 1E-14, unsigned count = 10000U) const; // conjugate gradients
	friend std::istream& operator>>(std::istream&, CRSMatrix&); // construct matrix from file
	friend std::ostream& operator<<(std::ostream&, CRSMatrix const &); // print matrix
	friend std::vector<REAL> operator*(CRSMatrix const &, std::vector<REAL> const &); // matrix-vector product
	friend CRSMatrix operator*(CRSMatrix const &, CRSMatrix const &); // matrix-matrix product
	friend std::vector<REAL> operator/(std::vector<REAL> const &, CRSMatrix const &); // solve A.x = b
};