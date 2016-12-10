#pragma once
#include "AbstractMultipliableMatrix.hpp"

template <typename T>
class DenseSquareMatrix : 
	public AbstractMultipliableMatrix<T> {
	T** _A;
	T*  _beg;
	// to be implemented
	T& _set(Index i, Index j) final { return _A[i][j]; }
	T  _get(Index i, Index j) const final { return _A[i][j]; };
public:
	explicit DenseSquareMatrix(Index n = 1);
	DenseSquareMatrix(std::initializer_list<std::initializer_list<T>> lst);
	DenseSquareMatrix(AbstractMatrix const &);
	~DenseSquareMatrix();
	DenseSquareMatrix& operator=(T const & val) final {
		for (Index i = 0; i < _w * _w; ++i) _beg[i] = val;
		return *this;
	}
	void mult(T const * by, T* result) const final;
	std::vector<T> GaussElimination(std::vector<T> const &);
};

// public methods

template <typename T>
DenseSquareMatrix<T>::DenseSquareMatrix(Index n) 
	: AbstractMatrix<T>(n, n)
	, _A(new double*[n])
	, _beg(new double[n * n]) {
	size_t i;
	for (i = 0; i < n; ++i) _A[i] = _beg + i * n;
	for (i = 0; i < n * n; ++i) _beg[i] = 0.;
}

template <typename T>
DenseSquareMatrix<T>::DenseSquareMatrix(std::initializer_list<std::initializer_list<T>> lst)
	: DenseSquareMatrix<T>(lst.size()) {
	Index i = 0, j = 0;
	for (auto const & row : lst) {
		for (auto const & val : row) _A[i][j++] = val;
		j = 0;
		++i;
	}
}

template <typename T>
DenseSquareMatrix<T>::DenseSquareMatrix(AbstractMatrix<T> const & A)
	: AbstractMatrix<T>(A.getOrder(), A.getOrder())
	, _A(new double*[A.getOrder()])
	, _beg(new double[A.getOrder() * A.getOrder()]) {
	auto n = A.getOrder();
	for (Index i = 0; i < n; ++i) {
		_A[i] = _beg + i * n;
		for (Index j = 0; j < n; ++j) _A[i][j] = A(i, j);
	}
}

template <typename T>
DenseSquareMatrix<T>::~DenseSquareMatrix() {
	delete[] _beg;
	delete[] _A;
}

template <typename T>
void DenseSquareMatrix<T>::mult(T const * by, T* result) const {
	for (Index i = 0; i < _h; ++i) {
		result[i] = 0.;
		for (Index j = 0; j < _w; ++j)
			result[i] += _A[i][j] * by[j];
	}
}

template <typename T>
std::vector<T> DenseSquareMatrix<T>::GaussElimination(std::vector<T> const & bConst) {
	auto n = getOrder();
	std::vector<T> x(n), b(bConst);
	int i;
	unsigned j, k, maxIndex;
	T* dummy;
	T max, mult;
	T sum;
	for (i = 0; i < n - 1; ++i) {
		// partial pivoting
		max = _A[i][i];
		maxIndex = i;
		for (j = i + 1; j < n; ++j)
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
		for (j = i + 1; j < n; ++j) {
			if (_A[j][i] == 0.) continue;
			mult = _A[j][i] / _A[i][i];
			// _A[j][i] = 0.;
			for (k = i + 1; k < n; ++k)
				_A[j][k] -= _A[i][k] * mult;
			b[j] -= b[i] * mult;
		}
	}
	for (i = n - 1; i >= 0; --i) {
		sum = 0.;
		for (j = i + 1; j < n; ++j)
			sum += _A[i][j] * x[j];
		x[i] = (b[i] - sum) / _A[i][i];
	}
	return x;
}