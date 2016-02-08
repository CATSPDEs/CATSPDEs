// declaration of CRSMatrix

#pragma once
#include "real.h"
#include <fstream>
#include <vector>
#include <cstddef> // size_t

class CRSMatrix { // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	size_t _n; // n := order of square matrix
	std::vector<size_t> _ptr, // ptr := locations in the val vector that start a row
						_col; // col := column indices of the elements in the val vector
	std::vector<REAL> _val; // val := values of the nonzero elements of the matrix (row-wised)
public:
	CRSMatrix(size_t, size_t); // order and number of nonzeros
	~CRSMatrix();
	size_t getOrder() const;
	void print() const;
	REAL operator()(size_t, size_t) const; // get value of element using indicies of traditional form of matrix
	REAL& operator()(size_t, size_t); // set —‘‘—
	friend std::istream& operator>>(std::istream&, CRSMatrix&); // construct from file
	friend std::vector<REAL> operator*(CRSMatrix const &, std::vector<REAL> const &); // matrix-vector product
	friend CRSMatrix operator*(CRSMatrix const &, CRSMatrix const &); // matrix-matrix product
	friend class SLAE;
};