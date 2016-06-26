#pragma once
#pragma once
#include "AbstractSparseMatrix.hpp"
#include "IRealMatrix.hpp"
#include "AdjacencyList.hpp"

class CSlRMatrix : public AbstractSparseMatrix<double>, public IRealMatrix {
	std::vector<double> _lval, // vector of elements of lower triangular part of matrix (raw by raw)
						_uval, // vector of elements of upper triangular part of matrix (raw by raw)
						_diag; // vector of diagonal elements (for FEM / FVM is always > 0, 
							   // so we store it explicitly)
	std::vector<size_t> _jptr, // vector of column (row) indicies and
						_iptr; // vector of raw indicies (see example to get what thes guys are)

	size_t _nnz() const { return _iptr[_n] + _n; } // look 2 strings above; “+ _n” because of _diag
	double& _set(size_t, size_t);
	double _get(size_t, size_t) const;
public:
	CSlRMatrix(size_t, size_t);
	~CSlRMatrix() {}
	std::vector<double> solve(std::vector<double> const &);
	std::vector<double> mult(std::vector<double> const &) const;
	std::vector<double> multt(std::vector<double> const &) const;
	std::istream& loadSparse(std::istream&); // look at implementation for istream / ostream structure
	std::ostream& saveSparse(std::ostream&) const;
	friend CSlRMatrix operator*(CSlRMatrix const &, CSlRMatrix const &);
	friend std::vector<double> operator&(CSlRMatrix const &, std::vector<double> const & );
};

