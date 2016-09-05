#pragma once
#include "AbstractSparseMatrix.hpp"
#include "AbstractMultipliableMatrix.hpp"
#include "IDecomposable.hpp"
#include "AdjacencyList.hpp"

template <typename T>
class SymmetricCSlRMatrix 
	: public AbstractSparseMatrix<T>
	, public AbstractMultipliableMatrix<T>
	, public IDecomposable<T> {
	// fancy name, yeah
	// stands for Symmetric Compressed Sparse (lower triangular) Row
	std::vector<double> _lval, // vector of elements of lower triangular part of matrix (raw by raw)
						_diag; // vector of diagonal elements (for FEM / FVM is always > 0, 
	                           // so we store it explicitly)
	std::vector<size_t> _jptr, // vector of column indicies and
	                    _iptr; // vector of row indicies (see example to get what thes guys are)
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
	// size of matrix is 6, so _iptr[6] = 8 is numb of elements in _lval or _jptr
	// makes sense!
	//
	// virtual methods to be implemented
	vector<T> _mult(vector<T> const &) const override;
	size_t _nnz() const override { return _iptr[_w] + _w; } // look 2 strings above; “+ _w” because of _diag
	T& _set(size_t, size_t) override;
	T  _get(size_t, size_t) const override;
	vector<T> _mult(vector<T> const &) const override;
public:
	SymmetricCSlRMatrix(size_t n, size_t nnz); // order of matrix and numb of nonzero elems
	SymmetricCSlRMatrix(AdjacencyList const &); // generate matrix portrait from adjacency list of mesh nodes 
	~SymmetricCSlRMatrix() {}
	// virtual methods to be implemented
	SymmetricCSlRMatrix& setZero() override;
	istream& loadSparse(istream&) override; // look at implementation for istream / ostream structure
	ostream& saveSparse(ostream&) const override;
	SymmetricCSlRMatrix& decompose() override;
	vector<T> forwardSubstitution (vector<T> const &) const override;
	vector<T> backwardSubstitution(vector<T> const &) const override;
};

