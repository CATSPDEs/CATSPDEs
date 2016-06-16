#pragma once
#include "AbstractSparseMatrix.hpp"
#include "IRealMatrix.hpp"

class CSRMatrix : public AbstractSparseMatrix<double>, public IRealMatrix { // matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	std::vector<size_t> _ptr, // ptr := locations in the val vector that start a row
						_col; // col := column indices of the elements in the val vector
	std::vector<double> _val; // val := values of the nonzero elements of the matrix (row-wised)
	// virtual methods to be implemented
	size_t _nnz() const { return _ptr[_n]; }
	double& _set(size_t, size_t);
	double _get(size_t, size_t) const;
public:
	CSRMatrix(size_t, size_t); // order and number of nonzeros
	~CSRMatrix() {}
	// virtual methods to be implemented
	std::vector<double> solve(std::vector<double> const &);
	std::vector<double> mult(std::vector<double> const &) const;
	CSRMatrix& setZero();
	std::istream& loadSparse(std::istream&);
	std::ostream& saveSparse(std::ostream&) const;
	// new methods
	std::vector<double> CG(std::vector<double> const &, double e = 1E-14, unsigned count = 10000U) const; // conjugate gradients
	friend CSRMatrix operator*(CSRMatrix const &, CSRMatrix const &); // matrix-matrix product
};