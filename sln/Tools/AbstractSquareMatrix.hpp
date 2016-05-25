#pragma once
#include <iostream>
#include <iomanip> // setprecision
#include <cstddef> // size_t
#include <stdexcept>
#include "vector.hpp" // vector operations (*, +, ...)

template <typename T>
class AbstractSquareMatrix { // abstract class for square matricies
	virtual T& _set(size_t, size_t) = 0; // set / get value of element using indicies of traditional form of matrix
	virtual T _get(size_t, size_t) const = 0;
protected: // interface for derivative classes
	size_t _n; // n := order of square matrix
public:
	AbstractSquareMatrix(size_t n) : _n(n) {
		if (_n < 1) throw std::out_of_range("order of matrix must be at least one");
	}
	virtual ~AbstractSquareMatrix() {}
	T operator()(size_t i, size_t j) const {
		if (i >= _n || j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _get(i, j);
	}
	T& operator()(size_t i, size_t j) { 
		if (i >= _n || j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _set(i, j); 
	}
	size_t getOrder() const { return _n; }
	std::ostream& save(std::ostream& output = std::cout) const {
		output << _n << '\n';
		output << std::setprecision(15) << std::scientific << std::showpos;
		for (size_t i = 0; i < _n; ++i) {
			for (size_t j = 0; j < _n; ++j)
				output << _get(i, j) << ' ';
			output << '\n';
		}
		return output;
	}
	std::istream& load(std::istream& input = std::cin) {
		T dummy;
		for (size_t i = 0; i < _n; ++i)
			for (size_t j = 0; j < _n; ++i)
				input >> dummy;
				try {
					_set(i, j) = dummy;
				} catch (std::invalid_argument const &) {} // you cannot set some values of sparse matrices because they are zeros
		return input;
	}
};