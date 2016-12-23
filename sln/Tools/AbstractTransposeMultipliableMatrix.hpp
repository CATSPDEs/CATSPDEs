#pragma once
#include "AbstractMultipliableMatrix.hpp"

template <typename T>
class AbstractTransposeMultipliableMatrix 
	: public AbstractMultipliableMatrix<T> {
	// enables v = A^T.u operation
	bool _t; 
public:
	AbstractTransposeMultipliableMatrix() : _t(false) {}
	AbstractTransposeMultipliableMatrix& t() { // flag for * method; see below
		_t = !_t;
		return *this;
	}
	virtual void multByTranspose(T const * by, T* result) const = 0;
	// tricky new version of simplified mult()
	// usage examples:
	// (1) ordinary multiplication
		// v = A * u;
	// (2) multiplication by transposed matrix (w/o actual transposing for sparse matrices, obviously)
		// v = A.t() * u;
	std::vector<T> operator*(std::vector<T> const & u) final {
		if (_t) {
			std::vector<T> v(_w, 0.);
			multByTranspose(u.data(), v.data());
			_t = false;
			return v;
		}
		std::vector<T> v(_h);
		mult(u.data(), v.data());
		return v;
	}
	// array version
	using AbstractMultipliableMatrix::operator*;
};