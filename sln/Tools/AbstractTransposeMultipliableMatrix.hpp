#pragma once
#include "AbstractMultipliableMatrix.hpp"

template <typename T>
class AbstractTransposeMultipliableMatrix 
	: public AbstractMultipliableMatrix<T> {
	// enables v = A^T.u operation
private:
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
	// (2) multiplication by transposed matrix (w/o actual transposing, obviously)
		// v = A.t() * u;
	vector<T> operator*(vector<T> const & u) final {
		if (_t) {
			vector<T> v(numbOfCols());
			multByTranspose(u.data(), v.data());
			_t = false;
			return v;
		}
		vector<T> v(numbOfRows());
		mult(u.data(), v.data());
		return v;
	}
};