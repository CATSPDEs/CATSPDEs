// declaration of BandMatrix

#pragma once
#include "AbstractSquareMatrix.hpp"

class BandMatrix : public AbstractSquareMatrix { // https://en.wikipedia.org/wiki/Band_matrix
	size_t _w, _p; // w := width of band (i.e. numb of non-zero diagonals), p := (w – 1) / 2
	bool _isDecomposed; // this flag shows whether matrix is LU decomposed 
	REAL* _A, ** _L, * _D, ** _U; // A = [ L | D | U ] := n × w band matrix (lower triangle, diagonal, upper triangle)
	bool _computeLU();
	void _forwardSubst(std::vector<REAL>&) const; // forward substitution; find vector y such that L.y = b
	void _backSubst(std::vector<REAL>&) const; // back substitution; find vector x such that U.x = y
public:
	BandMatrix(size_t, size_t); // create (allocate memory for) n × w matrix
	~BandMatrix(); // free memory
	// virtual methods to be implemented
	REAL& set(size_t, size_t);
	REAL get(size_t, size_t) const;
	std::vector<REAL> solve(std::vector<REAL> const &);
	std::istream& load(std::istream&);
	std::ostream& save(std::ostream&) const;
	std::vector<REAL> mult(std::vector<REAL> const &) const;
};