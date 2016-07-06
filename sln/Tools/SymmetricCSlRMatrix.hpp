#pragma once
#include "AbstractSparseMatrix.hpp"
#include "IRealMatrix.hpp"
#include "AdjacencyList.hpp"

class SymmetricCSlRMatrix : public AbstractSparseMatrix<double>, public IRealMatrix {
	// fancy name, yeah
	// stands for Symmetric Compressed Sparse (lower triangular) Row
	std::vector<double> _lval, // vector of elements of lower triangular part of matrix (raw by raw)
						_diag; // vector of diagonal elements (for FEM / FVM is always > 0, 
	                           // so we store it explicitly)
	std::vector<size_t> _jptr, // vector of column indicies and
	                    _iptr; // vector of raw indicies (see example to get what thes guys are)
	// example of 6 × 6 matrix:
	//		
	//  8 ─ ─ ─ ─ ┐
	//  × 1       |
	//  2 5 4	  |
	//  × × 3 2   |
	//  7 8 × 1 7 |
	//  × × 2 × 4 2
	//
	// (“×” denotes zero–element we do not store)
	//
	// _diag = { 8, 1, 4, 2, 7, 2 }
	// _lval = { 2, 5, 3, 7, 8, 1, 2, 4 }
	// _jptr = { 0, 1, 2, 0, 1, 3, 2, 4 }
	// _iptr = { 0, 0, 0, 2, 3, 6, 8 }
	// 
	// so, for one, 4th row (we count from zero) countains _iptr[5] – _iptr[4] = 6 – 3 = 3 elements, 
	// starting w/ element _lval[_iptr[4]] = _lval[3] = 7 in (_jptr[_iptr[4]] = _jptr[3] = 0)th column
	// size of matrix is 6, so _iptr[6 + 1] = 8 is numb of elements in _lval
	// makes sense!
	size_t _nnz() const { return _iptr[_n] + _n; } // look 2 strings above; “+ _n” because of _diag
	double& _set(size_t, size_t);
	double _get(size_t, size_t) const;
public:
	SymmetricCSlRMatrix(size_t, size_t);
	SymmetricCSlRMatrix(AdjacencyList const &); // generate matrix portrait from adjacency list of mesh nodes 
	~SymmetricCSlRMatrix() {}
	std::vector<double> solve(std::vector<double> const &);
	std::vector<double> mult(std::vector<double> const &) const;
	SymmetricCSlRMatrix& setZero();
	std::istream& loadSparse(std::istream&); // look at implementation for istream / ostream structure
	std::ostream& saveSparse(std::ostream&) const;
	SymmetricCSlRMatrix ILDL();
	std::vector<double> forwardSubst(std::vector<double> const &);
	friend SymmetricCSlRMatrix operator*(SymmetricCSlRMatrix const &, SymmetricCSlRMatrix const &);
};

