// implementation of CRSMatrix

#include <iostream>
#include "CRSMatrix.h"

CRSMatrix::CRSMatrix(size_t n, size_t nnz)
	: _n(n)
	, _ptr(n + 1)
	, _col(nnz)
	, _val(nnz) {
	_ptr[n] = nnz;
}

CRSMatrix::~CRSMatrix() {
}

// private methods

// public methods

size_t CRSMatrix::getOrder() const {
	return _n;
}

REAL CRSMatrix::operator()(size_t i, size_t j) const {
	for (size_t k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) break;
	return 0.;
}

REAL& CRSMatrix::operator()(size_t i, size_t j) {
	for (size_t k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) break;
	throw 0; // TODO "element w/ inicies (i, j) is zero and cannot be changed—you cannot change portrait of the matrix";
}

void CRSMatrix::print() const {
	size_t i, j;
	if (sizeof(REAL) == 4) std::cout.precision(6);
	else std::cout.precision(14);
	std::cout << std::scientific;
	for (i = 0; i < _n; ++i) {
		for (j = 0; j < _n; ++j) std::cout << (*this)(i, j) << ' ';
		std::cout << '\n';
	}
	std::cout << std::endl;
}

// friend functions

std::istream& operator>>(std::istream& input, CRSMatrix& A) {
	size_t i, max = A._n;
	for (i = 0; i < max; ++i)
		input >> A._ptr[i];
	max = A._ptr[max];
	for (i = 0; i < max; ++i)
		input >> A._col[i];
	for (i = 0; i < max; ++i)
		input >> A._val[i];
	return input;
}

std::vector<REAL> operator*(CRSMatrix const & A, std::vector<REAL> const & u) { // return product v = A.u
	std::vector<REAL> v(A._n, 0.);
	size_t i, j;
	for (i = 0; i < A._n; ++i)
		for (j = A._ptr[i]; j < A._ptr[i + 1]; ++j)
			v[i] += A._val[j] * u[A._col[j]];
	return v;
}

CRSMatrix operator*(CRSMatrix const & A, CRSMatrix const & B) { // return product N = A.B
	size_t n = A._n,
		   nnz = A._ptr[n];
	CRSMatrix N(n, nnz);
	// …
	return N;
}