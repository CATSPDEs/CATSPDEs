#pragma once
#include "CSlCMatrix.hpp"
#include "CSCMatrix.hpp"

/*
	Alexander Žilyakov, Jan 2017
*/

template <typename T>
class SaddlePointMatrix
	: public AbstractMultipliableMatrix<T> {
	T& _set(Index i, Index j) final {
		auto n = A11.getOrder();
		// A11 
		if (i < n && j < n) return A11(i, j);
		// A22
		if (n <= i && i < 2 * n &&
			n <= j && j < 2 * n) return A11(i - n, j - n);
		// to lower triangle
		if (i < j) std::swap(i, j);
		// not–A
		if (2 * n <= i) {
			// B11
			if (j < n) return B1(i - 2 * n, j);
			// B22
			if (j <= 2 * n) return B2(i - 2 * n, j - 2 * n);
		}
		throw std::invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
	}
	T  _get(Index i, Index j) const final {
		auto n = A11.getOrder();
		// A11 
		if (i < n && j < n) return A11(i, j);
		// A22
		if (n <= i && i < 2 * n &&
			n <= j && j < 2 * n) return A11(i - n, j - n);
		// to lower triangle
		if (i < j) std::swap(i, j);
		// not–A
		if (2 * n <= i) {
			// B11
			if (j < n) return B1(i - 2 * n, j);
			// B22
			if (j < 2 * n) return B2(i - 2 * n, j - n);
		}
		return 0.;
	}
public:
	//
	//  | A   B^T |
	//  | B   0   |
	//
	// A:
	//
	//  | A11   0   |
	//  | 0     A22 |
	//
	// A11 = A22
	//
	// B:
	//
	//  | B1  B2 |
	//

	CSlCMatrix<T> A11;
	CSCMatrix<double> B1, B2;
	SaddlePointMatrix(CSlCMatrix<T> const & a11, CSCMatrix<double> const & b1, CSCMatrix<double> const & b2) 
		: A11(a11), B1(b1), B2(b2), AbstractMatrix<T>(2 * a11.getOrder() + b1.numbOfRows())
	{}
	SaddlePointMatrix& operator=(T const & val) final {
		A11 = val; B1 = val; B2 = val;
		return *this;
	}
	void mult(T const * by, T* result) const final {
		// ...
	}
};