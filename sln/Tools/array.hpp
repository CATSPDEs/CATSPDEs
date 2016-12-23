#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip> // setprecision for i/o
#include "Index.hpp"

// sum two arrays
template <typename T, Index N>
std::array<T, N>& operator+=(std::array<T, N>& u, std::array<T, N> const & v) {
	std::transform(u.begin(), u.end(), v.begin(), u.begin(), std::plus<T>());
	return u;
}
template <typename T, Index N>
std::array<T, N> operator+(std::array<T, N> const & u, std::array<T, N> const & v) {
	std::array<T, N> t(u);
	return t += v;
}

// diff two arrays
template <typename T, Index N>
std::array<T, N>& operator-=(std::array<T, N>& u, std::array<T, N> const & v) {
	std::transform(u.begin(), u.end(), v.begin(), u.begin(), std::minus<T>());
	return u;
}
template <typename T, Index N>
std::array<T, N> operator-(std::array<T, N> const & u, std::array<T, N> const & v) {
	std::array<T, N> t(u);
	return t -= v;
}

// multiply by a scaler
template <typename T, Index N>
std::array<T, N>& operator*=(std::array<T, N>& u, T const & scaler) {
	std::for_each(u.begin(), u.end(), [&](T& elem) { elem *= scaler; });
	return u;
}
template <typename T, Index N>
std::array<T, N> operator*(T const & scaler, std::array<T, N> const & u) {
	std::array<T, N> t(u);
	return t *= scaler;
}

// divide by a scaler
template <typename T, Index N>
std::array<T, N>& operator/=(std::array<T, N>& u, T const & scaler) {
	return u *= 1 / scaler;
}
template <typename T, Index N>
std::array<T, N> operator/(std::array<T, N> const & u, T const & scaler) {
	std::array<T, N> t(u);
	return t /= scaler;
}

// i/o
template <typename T, Index N>
std::istream& operator>>(std::istream& input, std::array<T, N> & u) {
	std::for_each(u.begin(), u.end(), [&](T& elem) { input >> elem; });
	return input;
}
template <typename T, Index N>
std::ostream& operator<<(std::ostream& output, std::array<T, N> const & u) {
	output.precision(15/*4*/); // double precision
	output << std::scientific << std::showpos; // scientific notation (signed)
	std::for_each(u.begin(), u.end(), [&](T const & elem) { output << elem << ' '; });
	return output;
}

// dot product
template <typename T, Index N>
T operator*(std::array<T, N> const & u, std::array<T, N> const & v) {
	return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
}

// norm
template <typename T, Index N>
T norm(std::array<T, N> const & u) {
	return sqrt(u * u);
}

// middle point
template <typename T, Index N>
std::array<T, N> midNode(std::array<T, N> const & u, std::array<T, N> const & v) {
	return .5 * (u + v);
}

// זproduct (3D)
template <typename T>
inline std::array<T, 3> crossProduct(std::array<T, 3> const & u, std::array<T, 3> const & v) {
	return{
		u[1] * v[2] - u[2] * v[1],
		u[2] * v[0] - u[0] * v[2],
		u[0] * v[1] - u[1] * v[0]
	};
}

// זproduct (2D)
template <typename T>
inline std::array<T, 3> crossProduct(std::array<T, 2> const & u, std::array<T, 2> const & v) {
	return{ 0., 0, u[0] * v[1] - u[1] * v[0] };
}

// hash for arrays
// so one can use them as a key in maps or sets
namespace std {
	template <typename T, Index N>
	class hash<array<T, N>> {
	public:
		size_t operator()(array<T, N> const & u) const {
			return accumulate(
				next(u.begin()), 
				u.end(),
				u[0], // start with first element
				[](T const & a, T const & b) { return hash<T>()(a) ^ hash<T>()(b); }
			);
		}
	};
};

// we want std::array<T, 1> to be just T, not std::array<T, 1>
// so here goes some magic
template <typename T, Index D> struct ArrayTypedef       { typedef std::array<T, D> type; };
template <typename T>          struct ArrayTypedef<T, 1> { typedef T type; };
// and our smart array is
template <typename T, Index D>
using SmartArray = typename ArrayTypedef<T, D>::type;