#pragma once
#include <vector>
#include "Point.hpp"
#include "AbstractSquareMatrix.hpp"

typedef double (*fnctn)(Point const &); // here we are going to store... pointers to functions (defined on our points)! ~__~

class SymmetricMatrixOfFunctions : public AbstractSquareMatrix<fnctn> {
	std::vector<fnctn> _val;
public:
	explicit SymmetricMatrixOfFunctions(size_t n) : AbstractSquareMatrix(n), _val((n * n + n) / 2, nullptr) {}
	fnctn& set(size_t i, size_t j) {
		if (i > j) std::swap(i, j);
		if (j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
		return _val[j + _n * i - (i * i + i) / 2];
	}
};
