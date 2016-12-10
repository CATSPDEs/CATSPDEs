#pragma once
#include <fstream>
#include <iomanip> // setprecision
#include <stdexcept>
#include "SingletonLogger.hpp" // for logs
#include "Index.hpp"
#include "vector.hpp" // vector operations (*, +, ...)

template <typename T>
class AbstractMatrix { // abstract class for rect matricies
	virtual T& _set(Index, Index) = 0; // set / get value of element using indicies of traditional form of matrix
	virtual T  _get(Index, Index) const = 0;
protected: // interface for derivative classes
	Index _h, _w; // numb of rows and cols, _w × _h =: size of matrix
public:
	AbstractMatrix(Index h, Index w) 
		: _h(h)
		, _w(w) {
		if (_w < 1 || _h < 1) throw std::out_of_range("order of matrix must be at least one");
	}
	virtual ~AbstractMatrix() {}
	virtual AbstractMatrix& operator=(T const & val) = 0; // set all entries = @val (so we can “clear” matrix w/ A = 0.)
	T  operator()(Index i, Index j) const {
		if (i >= _h || j >= _w) throw std::out_of_range("matrix does not contain element w/ these indicies");
		return _get(i, j);
	}
	T& operator()(Index i, Index j) { 
		if (i >= _h || j >= _w) throw std::out_of_range("matrix does not contain element w/ these indicies");
		return _set(i, j); 
	}
	Index numbOfCols() const { return _w; }
	Index numbOfRows() const { return _h; }
	Index getOrder()   const { return _w; } // for square matrices (just for convenience)
	AbstractMatrix& import(std::istream& input = std::cin) {
		T dummy;
		setZero(); // clear matrix
		for (Index i = 0; i < _h; ++i)
			for (Index j = 0; j < _w; ++i) {
				input >> dummy;
				if (dummy) _set(i, j) = dummy; // to avoid exceptions if matrix portrait
				// doesn’t allow element w/ indicies (i, j)
			}
		return *this;
	}
	void export(std::ostream& output = std::cout) const {
		output << _h << ' ' << _w << '\n';
		output << std::setprecision(3/*15*/) << std::scientific << std::showpos;
		for (Index i = 0; i < _h; ++i) {
			for (Index j = 0; j < _w; ++j)
				output << _get(i, j) << ' ';
			output << '\n';
		}
	}
	AbstractMatrix& import(std::string const & inputStr) {
		std::ifstream input(inputStr);
		return import(input);
	}
	void export(std::string const & outputStr) const {
		std::ofstream output(outputStr);
		export(output);
	}
	virtual AbstractMatrix& operator*=(T const & val) {
		for (Index i = 0; i < _h; ++i) 
			for (Index j = 0; j < _w; ++j)
				_set(i, j) = val * _get(i, j);
		return *this;
	}
	AbstractMatrix& operator/=(T const & val) {
		return operator*=(1. / val);
	}
};