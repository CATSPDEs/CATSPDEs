// implementation of CRSMatrix

#include <numeric> // inner product
#include <cmath> // sqrt
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

std::vector<REAL> CRSMatrix::CG(std::vector<REAL> const & f, REAL e, unsigned iCount) const {
	unsigned n = 0;
	CRSMatrix const & A = *this; // synonim
	std::vector<REAL> x(_n, .0);
	std::vector<REAL> r(_n);
	std::vector<REAL> z;
	REAL a, b, dotProdR, dotProdF;
	r = f - A * x;
	z = r;
	dotProdR = r * r;
	dotProdF = f * f;
	while ((sqrt(dotProdR) / sqrt(dotProdF)) > e && n < iCount)
	{
		a = dotProdR / ((A*z) * z);
		x = x + a*z;
		r = r - a*(A*z);
		b = dotProdR;
		dotProdR = r * r;
		b = dotProdR / b;
		z = r + b*z;
		n++;
	}
	return x;
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

std::ostream& operator<<(std::ostream& output, CRSMatrix const & A) {
	size_t i, j, n = A._n;
	if (sizeof(REAL) == 4) output.precision(6);
	else output.precision(14);
	output << std::scientific;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) output << A(i, j) << ' ';
		output << '\n';
	}
	output << std::endl;
	return output;
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

std::vector<REAL> operator/(std::vector<REAL> const & b, CRSMatrix const & A) { // solve A.x = b
	// magic goes here
	// e.g., CG
	return A.CG(b);
}