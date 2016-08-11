#include <complex>
#include "CSRMatrix.hpp"

template <typename T>
CSRMatrix::CSRMatrix(size_t w, size_t h, size_t nnz)
	: AbstractSparseMatrix(w, h)
	, _iptr(h + 1)
	, _jptr(nnz)
	, _mval(nnz) {
	_iptr[h] = nnz;
}

// private methods

template <typename T>
T& CSRMatrix::_set(size_t i, size_t j) {
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _mval[k];
		else if (_jptr[k] > j) throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

template <typename T>
T CSRMatrix::_get(size_t i, size_t j) const {
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _mval[k];
		else if (_jptr[k] > j) return 0.;
	return 0.;
}

template <typename T>
vector<T> CSRMatrix::_mult(vector<T> const & u) const { // return product v = A.u
	vector<T> v(_h, 0.);
	size_t i, j;
	for (i = 0; i < _h; ++i)
		for (j = _iptr[i]; j < _iptr[i + 1]; ++j)
			v[i] += _mval[j] * u[_jptr[j]];
	return v;
}

template <typename T>
vector<T> CSRMatrix::_multByTranspose(vector<T> const & u) const { // return product v = A^T.u
	vector<T> v(_w, 0.);
	size_t i, j;
	for (i = 0; i < _h; ++i)
		for (j = _iptr[i]; j < _iptr[i + 1]; ++j)
			v[j] += _mval[j] * u[i];
	return v;
}

// public methods

template <typename T>
CSRMatrix& CSRMatrix::setZero() {
	fill(_mval.begin(), _mval.end(), 0.);
	return *this;
}

template <typename T>
istream& CSRMatrix::loadSparse(istream& input) {
	size_t i, max = _h;
	for (i = 0; i < max; ++i)
		input >> _iptr[i];
	max = _iptr[max];
	for (i = 0; i < max; ++i)
		input >> _jptr[i];
	for (i = 0; i < max; ++i)
		input >> _mval[i];
	return input;
}

template <typename T>
ostream& CSRMatrix::saveSparse(ostream& output) const {
	return output << _w << ' ' << _h << ' ' << _iptr[_h] << '\n'
				  << _iptr << '\n'
				  << _jptr << '\n'
				  << _mval;
}

// explicit instantiation

template class CSRMatrix<double>;
template class CSRMatrix<complex<double>>;
