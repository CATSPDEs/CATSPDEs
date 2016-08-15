#pragma once
#include "AbstractSparseMatrix.hpp"
#include "AbstractTransposeMultipliableMatrix.hpp"

template <typename T>
class CSRMatrix 
	: public AbstractSparseMatrix<T>
	, public AbstractTransposeMultipliableMatrix<T> { 
	// matrix is stored in comressed row storage (CSR for compressed sparse row) format (http://netlib.org/linalg/html_templates/node91.html)
	vector<T>      _mval; // values of the nonzero elements of the matrix (row-wised)
	vector<size_t> _iptr, // vector of column indicies and
	               _jptr; // vector of row indicies (see example to get what thes guys are)
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
	size_t _nnz() const final { return _iptr[_h]; }
	T& _set(size_t, size_t) final;
	T  _get(size_t, size_t) const final;
public:
	CSRMatrix(size_t h, size_t w, size_t nnz); // order and number of nonzeros
	// virtual methods to be implemented
	CSRMatrix& operator=(T const & val) final;
	void mult           (T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	istream& loadSparse(istream& from = cin) final;
	ostream& saveSparse(ostream& to   = cout) const final;
};

// implementation

template <typename T>
CSRMatrix<T>::CSRMatrix(size_t h, size_t w, size_t nnz)
	: AbstractMatrix<T>(h, w)
	, _iptr(h + 1)
	, _jptr(nnz)
	, _mval(nnz) {
	_iptr[h] = nnz;
}

// private methods

template <typename T>
T& CSRMatrix<T>::_set(size_t i, size_t j) {
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _mval[k];
		else if (_jptr[k] > j) throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

template <typename T>
T CSRMatrix<T>::_get(size_t i, size_t j) const {
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _mval[k];
		else if (_jptr[k] > j) return 0.;
	return 0.;
}

// public methods

template <typename T>
CSRMatrix<T>& CSRMatrix<T>::operator=(T const & val) {
	fill(_mval.begin(), _mval.end(), val);
	return *this;
}

template <typename T>
void CSRMatrix<T>::mult(T const * by, T* result) const {
	size_t i, j;
	for (i = 0; i < _h; ++i) {
		result[i] = 0.;
		for (j = _iptr[i]; j < _iptr[i + 1]; ++j)
			result[i] += _mval[j] * by[_jptr[j]];
	}
}

template <typename T>
void CSRMatrix<T>::multByTranspose(T const * by, T* result) const {
	size_t i, j;
	fill(result, result + numbOfCols(), 0.); // clear resulting vector
	for (i = 0; i < _h; ++i)
		for (j = _iptr[i]; j < _iptr[i + 1]; ++j)
			result[_jptr[j]] += _mval[j] * by[i];
}

template <typename T>
istream& CSRMatrix<T>::loadSparse(istream& input) {
	// stdin structure:
	// (1) _h + 1    elements of _iptr,
	// (2) _iptr[_h] elements of _jptr, and
	// (3) _iptr[_h] elements of _mval
	return input >> _iptr >> _jptr >> _mval;
}

template <typename T>
ostream& CSRMatrix<T>::saveSparse(ostream& output) const {
	return output << _h << ' ' << _w << ' ' << _iptr[_h] << '\n'
	              << _iptr << '\n'
		          << _jptr << '\n'
		          << _mval;
}
