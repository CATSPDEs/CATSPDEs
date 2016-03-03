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
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2)
{
	std::vector<T> res(v1.size());
	for (size_t i = 0;i < v1.size();i++)
		res[i] = v1[i] - v2[i];
	return res;
}

template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<T> const & u) {
	for (size_t i = 0; i < u.size(); ++i)
		output << u[i] << ' ';
	return output << std::endl;
}