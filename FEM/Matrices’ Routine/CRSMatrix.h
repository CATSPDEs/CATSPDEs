#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include "real.h"
#include <vector>

class CRSMatrix // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
{
	unsigned _n; // n := order of square matrix
	std::vector<REAL> _val; // val := values of the nonzero elements of the matrix (row-wised)
	std::vector<unsigned> _col, // col := column indices of the elements in the val vector
						  _ptr; // ptr := locations in the val vector that start a row
public:
	friend std::vector<REAL> operator*(const CRSMatrix&, const std::vector<REAL>&); // multiply matrix by vector
};

#endif