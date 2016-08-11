#pragma once
#include "AbstractSparseMatrix.hpp"
#include "ITransposeMultipliable.hpp"

template <typename T>
class CSRMatrix : public AbstractSparseMatrix<T>, public ITransposeMultipliable<T> { 
	// matrix is stored in comressed row storage format (http://netlib.org/linalg/html_templates/node91.html)
	vector<T>      _mval; // values of the nonzero elements of the matrix (row-wised)
	vector<size_t> _iptr, // vector of column indicies and
	               _jcol; // vector of row indicies (see example to get what thes guys are)
	// example of 3 × 6 matrix:
	//		
	//  5 × × 1 9 ×
	//  × 8 × × × 3
	//  1 4 × × 2 ×
	//
	// (“×” denotes zero–element we do not store)
	//
	// _mval = { 5, 1, 9, 8, 3, 1, 4, 2 }
	// _jptr = { 0, 3, 4, 1, 5, 0, 1, 4 }
	// _iptr = { 0, 3, 5, 8 }
	// 
	// so, for one, 1st row (we count from zero) countains _iptr[2] – _iptr[1] = 5 – 3 = 2 elements, 
	// starting w/ element _mval[_iptr[1]] = _mval[3] = 8 in (_jptr[_iptr[1]] = _jptr[3] = 1)st column
	// numb of rows is 3, so _iptr[3 + 1] = 8 is numb of elements in _mval and _jptr
	// makes sense!
	//
	// virtual methods to be implemented
	vector<T> _mult           (vector<T> const &) const override;
	vector<T> _multByTranspose(vector<T> const &) const override;
	size_t _nnz() const override { return _iptr[_h]; }
	T& _set(size_t, size_t) override;
	T  _get(size_t, size_t) const override;
public:
	CSRMatrix(size_t w, size_t h, size_t nnz); // order and number of nonzeros
	~CSRMatrix() {}
	// virtual methods to be implemented
	CSRMatrix& setZero() override;
	istream& loadSparse(istream&) override;
	ostream& saveSparse(ostream&) const override;
	// new methods
	//vector<double> CG(vector<double> const &, double e = 1E-14, unsigned count = 10000U) const; // conjugate gradients
	//friend CSRMatrix operator*(CSRMatrix const &, CSRMatrix const &); // matrix-matrix product
};