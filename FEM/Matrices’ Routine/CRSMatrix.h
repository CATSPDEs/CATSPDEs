// declaration of CRSMatrix

#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include "real.h"
#include <fstream>
#include <vector>

class CRSMatrix // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
{
	unsigned _n; // n := order of square matrix
	std::vector<REAL> _val; // val := values of the nonzero elements of the matrix (row-wised)
	std::vector<unsigned> _col, // col := column indices of the elements in the val vector
		_ptr; // ptr := locations in the val vector that start a row
public:
	CRSMatrix(unsigned);
	CRSMatrix(std::istream&); // construct from file 
	~CRSMatrix();
	unsigned getOrder() const;
	void print() const;
	REAL operator()(const unsigned, const unsigned) const; // get element using indicies of traditional form of matrix
	friend std::vector<REAL> operator*(const CRSMatrix&, const std::vector<REAL>&); // matrix-vector product
	friend CRSMatrix operator*(const CRSMatrix&, const CRSMatrix&); // matrix-matrix product
	friend class SLAE;
};

#endif