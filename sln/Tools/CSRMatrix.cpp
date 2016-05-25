// implementation of CSRMatrix

#include <cmath> // sqrt
#include "CSRMatrix.hpp"

CSRMatrix::CSRMatrix(size_t n, size_t nnz)
	: AbstractSparseMatrix(n)
	, _ptr(n + 1)
	, _col(nnz)
	, _val(nnz) {
	_ptr[n] = nnz;
}

// private methods

double& CSRMatrix::_set(size_t i, size_t j) {
	for (size_t k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) throw std::invalid_argument("portrait of sparse matrix doesn�t allow to change element w/ these indicies");
}

double CSRMatrix::_get(size_t i, size_t j) const {
	for (size_t k = _ptr[i]; k < _ptr[i + 1]; ++k)
		if (_col[k] == j) return _val[k];
		else if (_col[k] > j) return 0.;
}

// public methods

std::vector<double> CSRMatrix::mult(std::vector<double> const & u) const { // return product v = A.u
	std::vector<double> v(_n, 0.);
	size_t i, j;
	for (i = 0; i < _n; ++i)
		for (j = _ptr[i]; j < _ptr[i + 1]; ++j)
			v[i] += _val[j] * u[_col[j]];
	return v;
}

std::vector<double> CSRMatrix::solve(std::vector<double> const & b) {
	// magic goes here
	// e.g., CG
	return CG(b);
}

std::istream& CSRMatrix::loadSparse(std::istream& input) {
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

std::ostream& CSRMatrix::saveSparse(std::ostream& output) const {
	return output << _n << ' ' << _ptr[_n] << '\n'
				  << _ptr << '\n'
				  << _col << '\n'
				  << _val;
}

std::vector<double> CSRMatrix::CG(std::vector<double> const & f, double e, unsigned iCount) const {
	unsigned n = 0;
	CSRMatrix const & A = *this; // synonim
	std::vector<double> x(_n, .0);
	std::vector<double> r(_n);
	std::vector<double> z;
	double a, b, dotProdR, dotProdF;
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

// friends

CSRMatrix operator*(CSRMatrix const & A, CSRMatrix const & B) { // return product N = A.B
	size_t n = A._n,
		   nnz = A._ptr[n];
	CSRMatrix N(n, nnz);
	// �
	return N;
}