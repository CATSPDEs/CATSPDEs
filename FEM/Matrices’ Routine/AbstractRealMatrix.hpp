#pragma once
#include <iostream>
#include "vector.hpp"
#include "AbstractSquareMatrix.hpp"

class AbstractRealMatrix : public AbstractSquareMatrix<double> {
public:
	AbstractRealMatrix(size_t n) : AbstractSquareMatrix(n) {}
	virtual double get(size_t, size_t) const = 0;
	virtual std::vector<double> solve(std::vector<double> const &) = 0; // solve A.x = b
	virtual std::istream& load(std::istream&) = 0; // construct matrix from stdin (sparse form)
	virtual std::ostream& save(std::ostream&) const = 0; // save matrix to stdout (sparse form)
	virtual std::vector<double> mult(std::vector<double> const &) const = 0; // matrix-vector product

	void print(); // print matrix (non-sparse form)
};

// useful stuff 
std::istream& operator>>(std::istream&, AbstractRealMatrix&);
std::ostream& operator<<(std::ostream&, AbstractRealMatrix const &);
std::vector<double> operator*(AbstractRealMatrix const &, std::vector<double> const &);
std::vector<double> operator/(std::vector<double> const &, AbstractRealMatrix &);