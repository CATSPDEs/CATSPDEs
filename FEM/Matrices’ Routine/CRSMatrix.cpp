// implementation of CRSMatrix

#include <cmath> // sqrt
#include "CRSMatrix.hpp"

CRSMatrix::CRSMatrix(size_t n, size_t nnz)
	: AbstractSquareMatrix(n)
	, _ptr(n + 1)
	, _col(nnz)
	, _val(nnz) {
	_ptr[n] = nnz;
}

// public methods

REAL& CRSMatrix::set(size_t i, size_t j) {
	if (i >= _n || j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
	for (size_t k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

REAL CRSMatrix::get(size_t i, size_t j) const {
	if (i >= _n || j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
	for (size_t k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) return 0.;
}

std::vector<REAL> CRSMatrix::solve(std::vector<REAL> const & b) { // solve A.x = b
	// magic goes here
	// e.g., CG
	return CG(b);
}

std::istream& CRSMatrix::load(std::istream& input) {
	size_t i, max = _n;
	for (i = 0; i < max; ++i)
		input >> _ptr[i];
	max = _ptr[max];
	for (i = 0; i < max; ++i)
		input >> _col[i];
	for (i = 0; i < max; ++i)
		input >> _val[i];
	return input;
}

std::ostream& CRSMatrix::save(std::ostream& output) const {
	return output << _n << ' ' << _ptr[_n] << '\n'
				  << _ptr << '\n'
				  << _col << '\n'
				  << _val;
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

std::vector<REAL> CRSMatrix::mult(std::vector<REAL> const & u) const { // return product v = A.u
	std::vector<REAL> v(_n, 0.);
	size_t i, j;
	for (i = 0; i < _n; ++i)
		for (j = _ptr[i]; j < _ptr[i + 1]; ++j)
			v[i] += _val[j] * u[_col[j]];
	return v;
}

// friends

CRSMatrix operator*(CRSMatrix const & A, CRSMatrix const & B) { // return product N = A.B
	size_t n = A._n,
		   nnz = A._ptr[n];
	CRSMatrix N(n, nnz);
	// …
	return N;
}