#pragma once
#include <array>

using namespace std;

template <typename T, size_t N>
array<T, N>& operator*=(array<T, N>& u, T scaler) { // multiply by a scaler
	for (size_t i = 0; i < N; ++i)
		u[i] *= scaler;
	return u;
}

template <typename T, size_t N>
array<T, N> operator*(array<T, N> const & u, T scaler) {
	return array<T, N>(u) *= scaler;
}

template <typename T, size_t N>
array<T, N>& operator/=(array<T, N>& u, T scaler) { // divide by a scaler
	return u *= 1 / scaler;
}

template <typename T, size_t N>
array<T, N> operator/(array<T, N> const & u, T scaler) {
	return array<T, N>(u) /= scaler;
}