#pragma once
#include "AbstractMatrix.hpp"

template <typename T>
class AbstractMultipliableMatrix
	: public virtual AbstractMatrix<T> {
	// enables v = A.u operation
public:
	virtual void mult(T const * by, T* result) const = 0;
	// simplified version of mult() method
	virtual std::vector<T> operator*(std::vector<T> const & u) { // compute v = A.u
		std::vector<T> v(_h, 0.);
		mult(u.data(), v.data());
		return v;
	}
	// array version
	template <Index N>
	auto operator*(std::array<T, N> const & a) { 
		auto u = operator*(std::vector<T>(a.begin(), a.end()));
		std::array<T, N> v;
		std::copy_n(u.begin(), N, v.begin());
		return v;
	}
};