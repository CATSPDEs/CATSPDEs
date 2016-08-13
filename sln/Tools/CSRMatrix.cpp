#include <complex>
#include "CSRMatrix.hpp"

template <typename T>
CSRMatrix<T>::CSRMatrix(size_t w, size_t h, size_t nnz)
	: AbstractMatrix<T>(w, h)
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
	fill(result, result + numbOfRows(), 0.); // clear resulting vector
	for (i = 0; i < _h; ++i)
		for (j = _iptr[i]; j < _iptr[i + 1]; ++j)
			result[j] += _mval[j] * by[i];
}

template <typename T>
istream& CSRMatrix<T>::loadSparse(istream& input) {
	size_t i, max = _h;
	for (i = 0; i < max; ++i)
		input >> _iptr[i];
	max = _iptr[_h];
	for (i = 0; i < max; ++i)
		input >> _jptr[i];
	for (i = 0; i < max; ++i)
		input >> _mval[i];
	return input;
}

template <typename T>
ostream& CSRMatrix<T>::saveSparse(ostream& output) const {
	return output << _w << ' ' << _h << ' ' << _iptr[_h] << '\n'
				  << _iptr << '\n'
				  << _jptr << '\n'
				  << _mval;
}

// explicit instantiation

template class CSRMatrix<double>;
template class CSRMatrix<complex<double>>;
