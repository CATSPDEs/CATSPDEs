#pragma once
#include "AbstractSquareMatrix.hpp"

template <typename T>
class AbstractSparseMatrix : public AbstractSquareMatrix<T> {
	virtual size_t _nnz() const = 0; // number of nonzeros in sparse matrix
public:
	AbstractSparseMatrix(size_t n) : AbstractSquareMatrix(n) {}
	virtual std::istream& loadSparse(std::istream&) = 0; // construct matrix from stdin (sparse form)
	virtual std::ostream& saveSparse(std::ostream&) const = 0; // save matrix to stdout (sparse form)
	double sparsity() { return (double) _nnz() / _n; }
};

// useful stuff 
template <typename T>
std::istream& operator>>(std::istream& from, AbstractSparseMatrix<T>& A) { return A.loadSparse(from); }

template <typename T>
std::ostream& operator<<(std::ostream& to, AbstractSparseMatrix<T> const & A) { return A.saveSparse(to); }
