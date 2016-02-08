// declaration of CRSMatrix

#pragma once
#include "real.h"
#include <fstream>
#include <vector>
#include <cstddef> // size_t

class CRSMatrix { // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	size_t _n; // n := order of square matrix
	std::vector<REAL> _val; // val := values of the nonzero elements of the matrix (row-wised)
	std::vector<size_t> _col, // col := column indices of the elements in the val vector
						_ptr; // ptr := locations in the val vector that start a row
public:
	CRSMatrix(size_t);
	CRSMatrix(std::ifstream&); // construct from file 
	~CRSMatrix();
	size_t getOrder() const;
	void print() const;
	REAL operator()(size_t, size_t) const; // get value of element using indicies of traditional form of matrix
	REAL& operator()(size_t, size_t); // set —‘‘—
	friend std::vector<REAL> operator*(const CRSMatrix&, const std::vector<REAL>&); // matrix-vector product
	friend CRSMatrix operator*(const CRSMatrix&, const CRSMatrix&); // matrix-matrix product
	friend class SLAE;
};