#pragma once
#include <array>
#include <cmath>
#include <iostream>
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