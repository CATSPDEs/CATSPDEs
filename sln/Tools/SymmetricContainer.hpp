#pragma once
#include "AbstractSparseMatrix.hpp"

template <typename T>
class SymmetricContainer : public AbstractSparseMatrix<T> {
	std::vector<T> _val;
	// virtual methods to be implemented
	size_t _nnz() const { return (_n * _n + _n) / 2; }
	T& _set(size_t i, size_t j) {
		if (i > j) std::swap(i, j);
		return _val[j + _n * i - (i * i + i) / 2];
	}
	T _get(size_t i, size_t j) const {
		if (i > j) std::swap(i, j);
		return _val[j + _n * i - (i * i + i) / 2];
	}
public:
	explicit SymmetricContainer(size_t n) : AbstractSparseMatrix(n), _val(_nnz()) {}
	SymmetricContainer(size_t n, T val) : AbstractSparseMatrix(n), _val(_nnz(), val) {}
	// virtual methods to be implemented
	std::istream& loadSparse(std::istream& input) {
		return input >> _val;
	}
	std::ostream& saveSparse(std::ostream& output) const {
		return output << _n << '\n' << _val; 
	}
};