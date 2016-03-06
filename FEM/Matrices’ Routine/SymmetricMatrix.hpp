#pragma once
#include "AbstractRealMatrix.hpp"

class SymmetricMatrix : public AbstractRealMatrix {
	std::vector<double> _val;
public:
	explicit SymmetricMatrix(size_t n) : AbstractRealMatrix(n), _val((n * n + n) / 2, 0.) {}
	// to do: constructor from SymmetricMatrixOfFunctions at the given point
	double& set(size_t i, size_t j) {
		if (i > j) std::swap(i, j);
		if (j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _val[j + _n * i - (i * i + i) / 2];
	}
	double get(size_t i, size_t j) const {
		if (i > j) std::swap(i, j);
		if (j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _val[j + _n * i - (i * i + i) / 2];
	}
	// to do
	std::vector<double> solve(std::vector<double> const &) { return std::vector<double>(_n); }
	std::istream& load(std::istream& in) { return in; }
	std::ostream& save(std::ostream& output) const {
		return output << _n << '\n' << _val;
	}
	std::vector<double> mult(std::vector<double> const &) const { return std::vector<double>(_n); }
	SymmetricMatrix& identify() {
		std::fill(_val.begin(), _val.end(), 0.);
		for (size_t i = 0; i < _n; ++i) set(i, i) = 1.;
		return *this;
	}
};