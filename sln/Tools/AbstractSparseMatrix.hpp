#pragma once
#include "AbstractMatrix.hpp"

template <typename T>
class AbstractSparseMatrix 
	: public virtual AbstractMatrix<T> {
	// abstract class for sparse matrices
	virtual size_t _nnz() const = 0; // number of nonzeros in sparse matrix
public:
	virtual AbstractSparseMatrix& loadSparse(istream& from = cin) = 0; // construct matrix from stdin (sparse form)
	virtual void                  saveSparse(ostream& to   = cout) const = 0; // save matrix to stdout (sparse form)
	double sparsity() const { return (double)_nnz() / _w / _h; }
};

// useful stuff 
template <typename T>
istream& operator>>(istream& from, AbstractSparseMatrix<T>& A) { 
	A.loadSparse(from);
	return from; 
}

template <typename T>
ostream& operator<<(ostream& to, AbstractSparseMatrix<T> const & A) { 
	A.saveSparse(to);
	return to;
}
