#pragma once
#include <cstddef> // size_t
#include <stdexcept>

template <typename T>
class AbstractSquareMatrix { // abstract class for square matricies
protected: // interface for derivative classes
	size_t _n; // n := order of square matrix
	size_t _diff(size_t i, size_t j) const { return i > j ? i - j : j - i; }
public:
	AbstractSquareMatrix(size_t n) : _n(n) {
		if (_n < 1) throw std::out_of_range("order of matrix must be at least one");
	}
	virtual ~AbstractSquareMatrix() {}
	virtual T& set(size_t, size_t) = 0; // set / get value of element using indicies of traditional form of matrix
	T& operator()(size_t i, size_t j) { return set(i, j); }
	size_t getOrder() const { return _n; }
};