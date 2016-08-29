#pragma once
#include <iostream>
#include <iomanip> // setprecision
#include <cstddef> // size_t
#include <stdexcept>
#include "SingletonLogger.hpp" // for logs
#include "vector.hpp" // vector operations (*, +, ...)

template <typename T>
class AbstractMatrix { // abstract class for rect matricies
	virtual T& _set(size_t, size_t) = 0; // set / get value of element using indicies of traditional form of matrix
	virtual T  _get(size_t, size_t) const = 0;
protected: // interface for derivative classes
	size_t _h, _w; // numb of rows and cols, _w × _h =: size of matrix
public:
	AbstractMatrix(size_t h, size_t w) 
		: _h(h)
		, _w(w) {
		if (_w < 1 || _h < 1) throw out_of_range("order of matrix must be at least one");
	}
	virtual ~AbstractMatrix() {}
	virtual AbstractMatrix& operator=(T const & val) = 0; // set all entries = @val (so we can “clear” matrix w/ A = 0.)
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
	void save(ostream& output = cout) {
		output << _w << ' ' << _h << '\n';
		output << setprecision(6/*15*/) << scientific << showpos;
		for (size_t i = 0; i < _h; ++i) {
			for (size_t j = 0; j < _w; ++j)
				output << _get(i, j) << ' ';
			output << '\n';
		}
	}
	void load(istream& input = cin) {
		T dummy;
		setZero(); // clear matrix
		for (size_t i = 0; i < _h; ++i)
			for (size_t j = 0; j < _w; ++i) {
				input >> dummy;
				if (dummy) _set(i, j) = dummy; // to avoid exceptions if matrix portrait
				// doesn’t allow element w/ indicies (i, j)
			}
	}
};