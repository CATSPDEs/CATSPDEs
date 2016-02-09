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