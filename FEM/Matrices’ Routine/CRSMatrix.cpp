// implementation of CRSMatrix

#include <iostream>
#include "CRSMatrix.h"

CRSMatrix::CRSMatrix(unsigned n) :
	_n(n) {
	_ptr.resize(_n + 1);
}

CRSMatrix::CRSMatrix(std::istream &input) {
	unsigned i, m;
	input >> _n;
	_ptr.resize(_n + 1);
	for (i = 0; i <= _n; ++i)
		input >> _ptr[i];
	m = _ptr[_n];
	_val.resize(m);
	_col.resize(m);
	for (i = 0; i < m; ++i)
		input >> _val[i];
	for (i = 0; i < m; ++i)
		input >> _col[i];
}

CRSMatrix::~CRSMatrix() {
}

// private methods

// public methods

unsigned CRSMatrix::getOrder() const {
	return _n;
}

REAL CRSMatrix::operator()(const unsigned i, const unsigned j) const {
	unsigned k;
	for (k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) break;
	return 0.;
}

void CRSMatrix::print() const {
	unsigned i, j;
	if (sizeof(REAL) == 4) std::cout.precision(6);
	else std::cout.precision(14);
	std::cout << std::scientific;
	for (i = 0; i < _n; ++i) {
		for (j = 0; j < _n; ++j) std::cout << (*this)(i, j) << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// friend functions

std::vector<REAL> operator*(const CRSMatrix& A, const std::vector<REAL>& u) { // return product v = A.u
	std::vector<REAL> v(A._n, 0.);
	unsigned i, j;
	for (i = 0; i < A._n; ++i)
		for (j = A._ptr[i]; j < A._ptr[i + 1]; ++j)
			v[i] += A._val[j] * u[A._col[j]];
	return v;
}

CRSMatrix operator*(const CRSMatrix& A, const CRSMatrix& B) { // return product C = A.B
	unsigned n = A._n;
	CRSMatrix C(n);
	return C;
}