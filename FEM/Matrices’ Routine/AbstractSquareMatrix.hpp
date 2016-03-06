#pragma once
#include <cstddef> // size_t
#include <iostream>
#include <stdexcept>
#include "real.hpp"
#include "vector.hpp"

class AbstractSquareMatrix { // abstract class for square matricies
protected: // interface for derivative classes
	size_t _n; // n := order of square matrix
	size_t _diff(size_t i, size_t j) const { return i > j ? i - j : j - i; }
public:
	AbstractSquareMatrix(size_t n);
	virtual ~AbstractSquareMatrix() {}
	// pure virtual methods
	virtual REAL& set(size_t, size_t) = 0; // set / get value of element using indicies of traditional form of matrix
	virtual REAL get(size_t, size_t) const = 0;
	virtual std::vector<REAL> solve(std::vector<REAL> const &) = 0; // solve A.x = b
	virtual std::istream& load(std::istream&) = 0; // construct matrix from stdin (sparse form)
	virtual std::ostream& save(std::ostream&) const = 0; // save matrix to stdout (sparse form)
	virtual std::vector<REAL> mult(std::vector<REAL> const &) const = 0; // matrix-vector product

	REAL& operator()(size_t i, size_t j) { return set(i, j); };
	size_t getOrder() const { return _n; }
	void print(); // print matrix (non-sparse form)
};

// useful stuff 
std::istream& operator>>(std::istream&, AbstractSquareMatrix&);
std::ostream& operator<<(std::ostream&, AbstractSquareMatrix const &);
std::vector<REAL> operator*(AbstractSquareMatrix const &, std::vector<REAL> const &);
std::vector<REAL> operator/(std::vector<REAL> const &, AbstractSquareMatrix &);