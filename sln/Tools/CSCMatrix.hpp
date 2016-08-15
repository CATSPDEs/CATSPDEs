#pragma once
#include <algorithm> // for_each
#include "AbstractSparseMatrix.hpp"
#include "AbstractTransposeMultipliableMatrix.hpp"
#include "hb_io.hpp"

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
	CSCMatrix(ifstream& iHB); // construct from Harwell–Boeing file
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

template <typename T>
CSCMatrix<T>::CSCMatrix(ifstream& iHB)
	: AbstractMatrix<T>(1, 1) {
	int *colptr = nullptr; // why not size_t, HB_IO, why?..
	int indcrd;
	char *indfmt = nullptr;
	char *key = nullptr;
	char *mxtype = nullptr;
	int ncol;
	int neltvl;
	int nnzero;
	int nrhs;
	int nrhsix;
	int nrow;
	int ptrcrd;
	char *ptrfmt = nullptr;
	int rhscrd;
	char *rhsfmt = nullptr;
	char *rhstyp = nullptr;
	int *rowind = nullptr;
	char *title = nullptr;
	int totcrd;
	int valcrd;
	char *valfmt = nullptr;
	double *values = nullptr; // T = complex<double>?
	// read header info
	hb_header_read(iHB, &title, &key, &totcrd, &ptrcrd, &indcrd,
	               &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt,
	               &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix);
	if (mxtype[1] != 'R')
		throw invalid_argument("rect matrix is expected; use CSlR (SymmetricCSlR) class for unsymmetric matrices w/ symmetric pattern (symmetric matrices), respectively");
	if (mxtype[2] != 'A')
		throw invalid_argument("assembled matrix is expected");
	colptr = new int[ncol + 1];
	rowind = new int[nnzero];
	// read structure
	hb_structure_read(iHB, ncol, mxtype, nnzero, neltvl,
	                  ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind);
	values = new double[nnzero];
	// read values
	hb_values_read(iHB, valcrd, mxtype, nnzero, neltvl, valfmt, values);
	// construct matrix
	_h = nrow;
	_w = ncol;
	_jptr = vector<size_t>(colptr, colptr + ncol + 1);
	for_each(_jptr.begin(), _jptr.end(), [](size_t& i) { --i; });
	_iptr = vector<size_t>(rowind, rowind + nnzero);
	for_each(_iptr.begin(), _iptr.end(), [](size_t& i) { --i; });
	_mval = vector<T>(values, values + nnzero); // what if T = complex<double>?
	// free
	delete[] colptr;
	delete[] rowind;
	delete[] values;
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
	fill(result, result + _h, 0.); // clear resulting vector
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
