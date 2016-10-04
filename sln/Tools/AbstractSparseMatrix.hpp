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
	
	// in dense format, we store N = _w × _h elements; in sparse format we store some % of N elements (+ some extra data, such as column pointers, coords, band width etc.)
	// density() shows this % (not taking into account these extra data) 
	double density()  const { return (double) _nnz() / _w / _h; }
	double sparsity() const { return 1 - density(); }
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
