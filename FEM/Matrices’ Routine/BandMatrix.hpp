#pragma once
#include "AbstractSparseMatrix.hpp"
#include "IRealMatrix.hpp"

class BandMatrix : public AbstractSparseMatrix<double>, public IRealMatrix { // https://en.wikipedia.org/wiki/Band_matrix
	size_t _w, _p; // w := width of band (i.e. numb of non-zero diagonals), p := (w – 1) / 2
	bool _isDecomposed; // this flag shows whether matrix is LU decomposed 
	double* _A, ** _L, * _D, ** _U; // A = [ L | D | U ] := n × w band matrix (lower triangle, diagonal, upper triangle)
	bool _computeLU();
	void _forwardSubst(std::vector<double>&) const; // forward substitution; find vector y such that L.y = b
	void _backSubst(std::vector<double>&) const; // back substitution; find vector x such that U.x = y
	size_t _diff(size_t i, size_t j) const { return i > j ? i - j : j - i; }
	// virtual methods to be implemented
	size_t _nnz() const { return _n * _w; }
	double& _set(size_t, size_t);
	double _get(size_t, size_t) const;
public:
	BandMatrix(size_t, size_t); // create (allocate memory for) n × w matrix
	~BandMatrix(); // free memory
	// virtual methods to be implemented
	std::vector<double> solve(std::vector<double> const &);
	std::vector<double> mult(std::vector<double> const &) const;
	std::istream& loadSparse(std::istream&);
	std::ostream& saveSparse(std::ostream&) const;
};