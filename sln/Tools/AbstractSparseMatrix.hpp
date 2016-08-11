#pragma once
#include "AbstractMatrix.hpp"

template <typename T>
class AbstractSparseMatrix : public AbstractMatrix<T> {
	virtual size_t _nnz() const = 0; // number of nonzeros in sparse matrix
public:
	AbstractSparseMatrix(size_t w,size_t h) : AbstractMatrix(w,h) {}
	virtual istream& loadSparse(istream&) = 0; // construct matrix from stdin (sparse form)
	virtual ostream& saveSparse(ostream&) const = 0; // save matrix to stdout (sparse form)
	double sparsity() { return (double) _nnz() / _w / _h; }
};

// useful stuff 
template <typename T>
istream& operator>>(istream& from, AbstractSparseMatrix<T>& A) { return A.loadSparse(from); }

template <typename T>
ostream& operator<<(ostream& to, AbstractSparseMatrix<T> const & A) { return A.saveSparse(to); }
