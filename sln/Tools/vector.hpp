#pragma once
#include <vector>
#include <fstream> // i/o
#include <iomanip> // setprecision for i/o
#include "Index.hpp"

// dot product
template <typename T>
T operator*(std::vector<T> const & u, std::vector<T> const & v) {
	return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
}
template <typename T, Index N>
T operator*(std::vector<T> const & u, std::array<T, N> const & v) {
	return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
}
template <typename T, Index N>
T operator*(std::array<T, N> const & u, std::vector<T> const & v) {
	return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
}

// linear combination
template <typename T>
std::vector<T> operator*(std::vector<T> const & weights, std::vector<std::vector<T>> const & vectors) {
	return std::inner_product(vectors.begin(), vectors.end(), weights.begin(), 0.);
}
template <typename T, Index N>
std::array<T, N> operator*(std::vector<T> const & weights, std::vector<std::array<T, N>> const & vectors) {
	std::array<T, N> res;
	std::fill(res.begin(), res.end(), 0.);
	for (Index i = 0; i < weights.size(); ++i)
		res += weights[i] * vectors[i];
	return res;
}

// sum two vectors
template <typename T>
std::vector<T>& operator+=(std::vector<T>& u, std::vector<T> const & v) {
	std::transform(u.begin(), u.end(), v.begin(), u.begin(), std::plus<T>());
	return u;
}
template <typename T>
std::vector<T> operator+(std::vector<T> const & u, std::vector<T> const & v) {
	auto t(u);
	return t += v;
}

// diff two vectors
template <typename T>
std::vector<T>& operator-=(std::vector<T>& u, std::vector<T> const & v) {
	std::transform(u.begin(), u.end(), v.begin(), u.begin(), std::minus<T>());
	return u;
}
template <typename T>
std::vector<T> operator-(std::vector<T> const & u, std::vector<T> const & v) {
	auto t(u);
	return t -= v;
}

// multiply by a scaler
template <typename T>
std::vector<T>& operator*=(std::vector<T>& u, T const & scaler) {
	std::for_each(u.begin(), u.end(), [&](T& elem) { elem *= scaler; });
	return u;
}
template <typename T>
std::vector<T> operator*(T const & scaler, std::vector<T> const & u) {
	auto t(u);
	return t *= scaler;
}

// i/o
template <typename T>
std::istream& operator>>(std::istream& input, std::vector<T> & u) {
	std::for_each(u.begin(), u.end(), [&](T& elem) { input >> elem; });
	return input;
}
template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<T> const & u) {
	output.precision(15); // double precision
	output << std::scientific << std::showpos; // scientific notation (signed)
	std::for_each(u.begin(), u.end(), [&](T const & elem) { output << elem << ' '; });
	return output;
}
template <typename T>
void import(std::vector<T>& u, std::string const & importString) {
	std::ifstream import(importString);
	import >> u;
}
template <typename T>
void export(std::vector<T> const & u, std::string const & outputString) {
	std::ofstream output(outputString);
	output << u;
}


// for std::vector of std::vectors
template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<std::vector<T>> const & u) {
	output.precision(15); // double precision
	output << std::scientific << std::showpos;
	for (Index i = 0; i < u.size(); ++i) {
		for (Index j = 0; j < u[i].size(); ++j)
			output << u[i][j] << ' ';
		output << '\n';
	}
	return output;
}

template <typename T, Index N>
std::ostream& operator<<(std::ostream& output, std::vector<std::array<T, N>> const & u) {
	output.precision(15); // double precision
	output << std::scientific << std::showpos;
	for (Index i = 0; i < u.size(); ++i) {
		for (Index j = 0; j < N; ++j)
			output << u[i][j] << ' ';
		output << '\n';
	}
	return output;
}

template <typename T>
T norm(std::vector<T> const & v) {
	return sqrt(v * v);
}