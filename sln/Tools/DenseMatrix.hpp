#pragma once
#include <algorithm> // for_each
#include "AbstractTransposeMultipliableMatrix.hpp"

template <typename T>
class DenseMatrix :
	public AbstractTransposeMultipliableMatrix<T> {
	std::vector<std::vector<T>> _A;
	// to be implemented
	T& _set(Index i, Index j) final { return _A[i][j]; }
	T  _get(Index i, Index j) const final { return _A[i][j]; };
public:
	// create h × w zero matrix
	explicit DenseMatrix(Index h, Index w)
		: AbstractMatrix<T>(h, w)
		, _A(h, std::vector<T>(w, 0.))
	{}
	// create matrix from ini lists, e.g. from { {11, 12}, {21, 22} }
	DenseMatrix(std::initializer_list<std::initializer_list<T>> const &);
	// densify any matrix
	DenseMatrix(AbstractMatrix const &);
	DenseMatrix& operator=(T const & val) final {
		std::for_each(_A.begin(), _A.end(), [&](auto& row) {
			std::fill(row.begin(), row.end(), val);
		});
		return *this;
	}
	void mult           (T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	// Gauss elimination w/ partial pivoting
	std::vector<T> GaussElimination(std::vector<T> const &);
};

// public methods

template <typename T>
DenseMatrix<T>::DenseMatrix(std::initializer_list<std::initializer_list<T>> const & iniList)
	: AbstractMatrix<T>(iniList.size(), (*iniList.begin()).size()) {
	for (auto const & row : iniList) {
		if (row.size() != _w) throw std::invalid_argument("ini list does not represent a matrix");
		_A.emplace_back(row);
	}
}

template <typename T>
DenseMatrix<T>::DenseMatrix(AbstractMatrix<T> const & A)
	: DenseMatrix<T>(A.numbOfRows(), A.numbOfCols()) {
	for (Index i = 0; i < _h; ++i) 
		for (Index j = 0; j < _w; ++j) 
			_A[i][j] = A(i, j);
}

template <typename T>
void DenseMatrix<T>::mult(T const * by, T* result) const {
	for (Index i = 0; i < _h; ++i) 
		for (Index j = 0; j < _w; ++j)
			result[i] += _A[i][j] * by[j];
}

template <typename T>
void DenseMatrix<T>::multByTranspose(T const * by, T* result) const {
	for (Index j = 0; j < _w; ++j)
		for (Index i = 0; i < _h; ++i)
			result[j] += _A[i][j] * by[i];
}

template <typename T>
std::vector<T> DenseMatrix<T>::GaussElimination(std::vector<T> const & bConst) {
	if (_w != _h) throw std::invalid_argument("matrix must be square");
	auto n = getOrder();
	std::vector<T> x(n), b(bConst);
	SignedIndex i;
	Index j, k, maxIndex;
	T maxElem, mult, sum;
	for (i = 0; i < n - 1; ++i) {
		// partial pivoting
		maxElem = _A[maxIndex = i][i];
		for (k = i + 1; k < n; ++k) // find row w / max element
			if (_A[k][i] > maxElem) 
				maxElem = _A[maxIndex = k][i];
		_A[i].swap(_A[maxIndex]);
		std::swap(b[i], b[maxIndex]);
		// here Gauss goes!
		for (k = i + 1; k < n; ++k) {
			if (_A[k][i] == 0.) continue;
			mult = _A[k][i] / _A[i][i];
			for (j = i + 1; j < n; ++j)
				_A[k][j] -= _A[i][j] * mult;
			b[k] -= b[i] * mult;
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

//template <typename T>
//class DenseSquareMatrix : 
//	public AbstractMultipliableMatrix<T> {
//	T** _A;
//	T*  _beg;
//	// to be implemented
//	T& _set(Index i, Index j) final { return _A[i][j]; }
//	T  _get(Index i, Index j) const final { return _A[i][j]; };
//public:
//	explicit DenseSquareMatrix(Index n = 1);
//	DenseSquareMatrix(std::initializer_list<std::initializer_list<T>> lst);
//	DenseSquareMatrix(AbstractMatrix const &);
//	~DenseSquareMatrix();
//	DenseSquareMatrix& operator=(T const & val) final {
//		for (Index i = 0; i < _w * _w; ++i) _beg[i] = val;
//		return *this;
//	}
//	void mult(T const * by, T* result) const final;
//	std::vector<T> GaussElimination(std::vector<T> const &);
//};
//
//// public methods
//
//template <typename T>
//DenseSquareMatrix<T>::DenseSquareMatrix(Index n) 
//	: AbstractMatrix<T>(n, n)
//	, _A(new double*[n])
//	, _beg(new double[n * n]) {
//	size_t i;
//	for (i = 0; i < n; ++i) _A[i] = _beg + i * n;
//	for (i = 0; i < n * n; ++i) _beg[i] = 0.;
//}
//
//template <typename T>
//DenseSquareMatrix<T>::DenseSquareMatrix(std::initializer_list<std::initializer_list<T>> lst)
//	: DenseSquareMatrix<T>(lst.size()) {
//	Index i = 0, j = 0;
//	for (auto const & row : lst) {
//		for (auto const & val : row) _A[i][j++] = val;
//		j = 0;
//		++i;
//	}
//}
//
//template <typename T>
//DenseSquareMatrix<T>::DenseSquareMatrix(AbstractMatrix<T> const & A)
//	: AbstractMatrix<T>(A.getOrder(), A.getOrder())
//	, _A(new double*[A.getOrder()])
//	, _beg(new double[A.getOrder() * A.getOrder()]) {
//	auto n = A.getOrder();
//	for (Index i = 0; i < n; ++i) {
//		_A[i] = _beg + i * n;
//		for (Index j = 0; j < n; ++j) _A[i][j] = A(i, j);
//	}
//}
//
//template <typename T>
//DenseSquareMatrix<T>::~DenseSquareMatrix() {
//	delete[] _beg;
//	delete[] _A;
//}
//
//template <typename T>
//void DenseSquareMatrix<T>::mult(T const * by, T* result) const {
//	for (Index i = 0; i < _h; ++i) {
//		result[i] = 0.;
//		for (Index j = 0; j < _w; ++j)
//			result[i] += _A[i][j] * by[j];
//	}
//}
//
//template <typename T>
//std::vector<T> DenseSquareMatrix<T>::GaussElimination(std::vector<T> const & bConst) {
//	auto n = getOrder();
//	std::vector<T> x(n), b(bConst);
//	int i;
//	unsigned j, k, maxIndex;
//	T* dummy;
//	T max, mult;
//	T sum;
//	for (i = 0; i < n - 1; ++i) {
//		// partial pivoting
//		max = _A[i][i];
//		maxIndex = i;
//		for (j = i + 1; j < n; ++j)
//			if (_A[j][i] > max) {
//				max = _A[j][i];
//				maxIndex = j;
//			}
//		if (maxIndex != i) {
//			dummy = _A[i];
//			_A[i] = _A[maxIndex];
//			_A[maxIndex] = dummy;
//			max = b[i]; // max is also “dummy” here
//			b[i] = b[maxIndex];
//			b[maxIndex] = max;
//		}
//		// here Gauss goes!
//		for (j = i + 1; j < n; ++j) {
//			if (_A[j][i] == 0.) continue;
//			mult = _A[j][i] / _A[i][i];
//			// _A[j][i] = 0.;
//			for (k = i + 1; k < n; ++k)
//				_A[j][k] -= _A[i][k] * mult;
//			b[j] -= b[i] * mult;
//		}
//	}
//	for (i = n - 1; i >= 0; --i) {
//		sum = 0.;
//		for (j = i + 1; j < n; ++j)
//			sum += _A[i][j] * x[j];
//		x[i] = (b[i] - sum) / _A[i][i];
//	}
//	return x;
//}