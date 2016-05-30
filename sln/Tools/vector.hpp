#pragma once
#include <vector>

template <typename T>
std::vector<T> operator*(T c, std::vector<T> const & v)
{
	std::vector<T> res(v.size());
	for (size_t i = 0; i < v.size(); i++)
		res[i] = c*v[i];
	return res;
}

template <typename T>
T operator*(std::vector<T> const & u, std::vector<T> const & v) { // dot product
	T dotProduct = 0;
	for (size_t i = 0; i < u.size(); ++i)
		dotProduct += u[i] * v[i];
	return dotProduct;
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2)
{
	std::vector<T> res(v1.size());
	for (size_t i = 0;i < v1.size();i++)
		res[i] = v1[i] + v2[i];
	return res;
}

template <typename T>
std::vector<T> operator-(std::vector<T> const & v1, std::vector<T> const & v2) {
	std::vector<T> res(v1.size());
	for (size_t i = 0; i < v1.size(); i++)
		res[i] = v1[i] - v2[i];
	return res;
}

template <typename T>
std::istream& operator>>(std::istream& input, std::vector<T> & u) {
	for (size_t i = 0; i < u.size(); ++i)
		input >> u[i];
	return input;
}

template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<T> const & u) {
	output.precision(15); // double precision
	output << std::scientific << std::showpos;
	for (size_t i = 0; i < u.size(); ++i)
		output << u[i] << '\n';
	return output;
}

template <typename T>
T norm(std::vector<T> const & v)
{
	return sqrt(v*v);
}