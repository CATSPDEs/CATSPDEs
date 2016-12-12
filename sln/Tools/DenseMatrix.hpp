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