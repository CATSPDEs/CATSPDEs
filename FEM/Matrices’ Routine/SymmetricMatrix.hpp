#pragma once
#include "AbstractSquareMatrix.hpp"

class SymmetricMatrix : public AbstractSquareMatrix {
	std::vector<REAL> _val;
public:
	SymmetricMatrix(size_t n) : AbstractSquareMatrix(n), _val((n * n + n) / 2, 0.) {}
	REAL& set(size_t i, size_t j) {
		if (i > j) std::swap(i, j);
		if (j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _val[j + _n * i - (i * i + i) / 2];
	}
	REAL get(size_t i, size_t j) const {
		if (i > j) std::swap(i, j);
		if (j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _val[j + _n * i - (i * i + i) / 2];
	}
	// to do
	std::vector<REAL> solve(std::vector<REAL> const &) { return std::vector<REAL>(_n); }
	std::istream& load(std::istream& in) { return in; }
	std::ostream& save(std::ostream& output) const {
		return output << _n << '\n' << _val;
	}
	std::vector<REAL> mult(std::vector<REAL> const &) const { return std::vector<REAL>(_n); }
	SymmetricMatrix& identify() {
		std::fill(_val.begin(), _val.end(), 0.);
		//for (size_t i = 0; i < _n; ++i) set(i, i) = 1.;
		return *this;
	}
};