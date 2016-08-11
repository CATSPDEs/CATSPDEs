#pragma once
#include <iostream>
#include <iomanip> // setprecision
#include <cstddef> // size_t
#include <stdexcept>
#include "vector.hpp" // vector operations (*, +, ...)

template <typename T>
class AbstractMatrix { // abstract class for rect matricies
	virtual T& _set(size_t, size_t) = 0; // set / get value of element using indicies of traditional form of matrix
	virtual T  _get(size_t, size_t) const = 0;
protected: // interface for derivative classes
	size_t _w, _h; // numb of rows and cols, _w × _h =: size of matrix
public:
	explicit AbstractMatrix(size_t w, size_t h) 
		: _w(w)
		, _h(h) {
		if (_w < 1 || _h < 1) throw out_of_range("order of matrix must be at least one");
	}
	virtual ~AbstractMatrix() {}
	virtual AbstractMatrix& setZero() = 0; // set all entries = “0” (whatever “0” for T means)
	T operator()(size_t i, size_t j) const {
		if (i >= _h || j >= _w) throw out_of_range("matrix doesn’t contain element w/ these indicies");
		return _get(i, j);
	}
	T& operator()(size_t i, size_t j) { 
		if (i >= _h || j >= _w) throw out_of_range("matrix doesn’t contain element w/ these indicies");
		return _set(i, j); 
	}
	size_t numbOfCols() const { return _w; }
	size_t numbOfRows() const { return _h; }
	size_t getOrder()   const { return _w; } // for square matrices (just for convenience)
	ostream& save(ostream& output = cout) const {
		output << _n << '\n';
		output << setprecision(6/*15*/) << scientific << showpos;
		for (size_t i = 0; i < _h; ++i) {
			for (size_t j = 0; j < _w; ++j)
				output << _get(i, j) << ' ';
			output << '\n';
		}
		return output;
	}
	istream& load(istream& input = cin) {
		T dummy;
		for (size_t i = 0; i < _h; ++i)
			for (size_t j = 0; j < _w; ++i) {
				input >> dummy;
				if (dummy) _set(i, j) = dummy;
			}
		return input;
	}
};