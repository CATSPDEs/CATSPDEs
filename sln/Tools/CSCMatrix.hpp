#pragma once
#include "AbstractSparseMatrix.hpp"
#include "AbstractTransposeMultipliableMatrix.hpp"

template <typename T> class CSCMatrix;
// matrix is stored in CSC (compressed sparse column) format…

template <typename T> using HBMatrix = CSCMatrix<T>;
// …aka Harwell–Boeing format…

template <typename T>
class CSCMatrix
	: public AbstractSparseMatrix<T>
	, public AbstractTransposeMultipliableMatrix<T> {
	// …it’s similar to CSR (check out CSRMatrix.hpp for details) except 
	// _mval vector contains nozero matrix values col by col, not row by row
	// so we omit details here
	vector<T>      _mval; 
	vector<size_t> _jptr, 
	               _iptr; 
	// virtual methods to be implemented
	size_t _nnz() const final { return _jptr[_w]; }
	T& _set(size_t, size_t) final;
	T  _get(size_t, size_t) const final;
public:
	CSCMatrix(size_t h, size_t w, size_t nnz);
	// virtual methods to be implemented
	CSCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	istream& loadSparse(istream& from = cin) final;
	ostream& saveSparse(ostream& to   = cout) const final;
};

// implementation

template <typename T>
CSCMatrix<T>::CSCMatrix(size_t h, size_t w, size_t nnz)
	: AbstractMatrix<T>(h, w)
	, _jptr(w + 1)
	, _iptr(nnz)
	, _mval(nnz) {
	_jptr[w] = nnz;
}

// private methods

template <typename T>
T& CSCMatrix<T>::_set(size_t i, size_t j) {
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _mval[k];
		else if (_iptr[k] > i) throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

template <typename T>
T CSCMatrix<T>::_get(size_t i, size_t j) const {
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _mval[k];
		else if (_iptr[k] > i) return 0.;
	return 0.;
}

// public methods

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::operator=(T const & val) {
	fill(_mval.begin(), _mval.end(), val);
	return *this;
}

template <typename T>
void CSCMatrix<T>::mult(T const * by, T* result) const {
	size_t i, j;
	fill(result, result + numbOfRows(), 0.); // clear resulting vector
	for (j = 0; j < _w; ++j)
		for (i = _jptr[j]; i < _jptr[j + 1]; ++i)
			result[_iptr[i]] += _mval[i] * by[j];
}

template <typename T>
void CSCMatrix<T>::multByTranspose(T const * by, T* result) const {
	size_t i, j;
	for (j = 0; j < _w; ++j) {
		result[j] = 0.;
		for (i = _jptr[j]; i < _jptr[j + 1]; ++i)
			result[j] += _mval[i] * by[_iptr[i]];
	}
}

template <typename T>
istream& CSCMatrix<T>::loadSparse(istream& input) {
	// stdin structure:
	// (1) _w + 1    elements of _jptr,
	// (2) _jptr[_w] elements of _iptr, and
	// (3) _jptr[_w] elements of _mval
	return input >> _jptr >> _iptr >> _mval;
}

template <typename T>
ostream& CSCMatrix<T>::saveSparse(ostream& output) const {
	return output << _h << ' ' << _w << ' ' << _jptr[_w] << '\n'
	              << _jptr << '\n'
	              << _iptr << '\n'
	              << _mval;
}
