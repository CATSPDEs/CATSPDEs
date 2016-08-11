#pragma once
#include <vector>

using namespace std;

template <typename T>
vector<T> operator*(T c, vector<T> const & v)
{
	vector<T> res(v.size());
	for (size_t i = 0; i < v.size(); i++)
		res[i] = c*v[i];
	return res;
}

template <typename T>
T operator*(vector<T> const & u, vector<T> const & v) { // dot product
	T dotProduct = 0;
	for (size_t i = 0; i < u.size(); ++i)
		dotProduct += u[i] * v[i];
	return dotProduct;
}

template <typename T>
vector<T> operator+(const vector<T>& v1, const vector<T>& v2)
{
	vector<T> res(v1.size());
	for (size_t i = 0;i < v1.size();i++)
		res[i] = v1[i] + v2[i];
	return res;
}

template <typename T>
vector<T> operator-(vector<T> const & v1, vector<T> const & v2) {
	vector<T> res(v1.size());
	for (size_t i = 0; i < v1.size(); i++)
		res[i] = v1[i] - v2[i];
	return res;
}

template <typename T>
istream& operator>>(istream& input, vector<T> & u) {
	for (size_t i = 0; i < u.size(); ++i)
		input >> u[i];
	return input;
}

template <typename T>
ostream& operator<<(ostream& output, vector<T> const & u) {
	output.precision(15); // double precision
	output << scientific << showpos;
	for (size_t i = 0; i < u.size(); ++i)
		output << u[i] << '\n';
	return output;
}

// for vector of vectors
template <typename T>
ostream& operator<<(ostream& output, vector<vector<T>> const & u) {
	output.precision(15); // double precision
	output << scientific << showpos;
	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < u[i].size(); ++j)
			output << u[i][j] << ' ';
		output << '\n';
	}
	return output;
}

template <typename T, size_t N>
ostream& operator<<(ostream& output, vector<array<T, N>> const & u) {
	output.precision(15); // double precision
	output << scientific << showpos;
	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < N; ++j)
			output << u[i][j] << ' ';
		output << '\n';
	}
	return output;
}

template <typename T>
T norm(vector<T> const & v)
{
	return sqrt(v*v);
}