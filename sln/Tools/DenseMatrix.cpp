#include "DenseMatrix.hpp"

// private methods

double& DenseMatrix::_set(size_t i, size_t j) {
	return _A[i][j];
}

double DenseMatrix::_get(size_t i, size_t j) const {
	return _A[i][j];
}

// public methods

DenseMatrix::DenseMatrix(size_t n) :
	AbstractMatrix(n), _A(new double*[n]), _beg(new double[n * n]) {
	size_t i;
	for (i = 0; i < n; ++i)
		_A[i] = _beg + i * n;
	for (i = 0; i < n * n; ++i)
		_beg[i] = 0.;
}

DenseMatrix::~DenseMatrix() {
	if (_n) {
		delete[] _beg;
		delete[] _A;
	}
}

// void Matrix::solveGauss(double* x, double* b) { // x <— A.? = b
std::vector<double> DenseMatrix::solve(std::vector<double> const & bConst) {
	std::vector<double> x(_n), b(bConst);
	int i;
	unsigned j, k, maxIndex;
	double* dummy;
	double max, mult;
	double sum;
	for (i = 0; i < _n - 1; ++i) {
		// partial pivoting
		max = _A[i][i];
		maxIndex = i;
		for (j = i + 1; j < _n; ++j) 
			if (_A[j][i] > max) {
				max = _A[j][i];
				maxIndex = j;
			}
		if (maxIndex != i) {
			dummy = _A[i];
			_A[i] = _A[maxIndex];
			_A[maxIndex] = dummy;
			max = b[i]; // max is also “dummy” here
			b[i] = b[maxIndex];
			b[maxIndex] = max;
		}
		// here Gauss goes!
		for (j = i + 1; j < _n; ++j) {
			if (_A[j][i] == 0.) continue;
			mult = _A[j][i] / _A[i][i];
			// _A[j][i] = 0.;
			for (k = i + 1; k < _n; ++k) 
				_A[j][k] -= _A[i][k] * mult;
			b[j] -= b[i] * mult;
		}
	}
	for (i = _n - 1; i >= 0; --i) {
		sum = 0.;
		for (j = i + 1; j < _n; ++j)
			sum += _A[i][j] * x[j];
		x[i] = (b[i] - sum) / _A[i][i];
	}
	return x;
}

std::vector<double> DenseMatrix::mult(std::vector<double> const & v) const {
	return std::vector<double>(_n);
}