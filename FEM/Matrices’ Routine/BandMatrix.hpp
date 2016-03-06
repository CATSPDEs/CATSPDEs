// declaration of BandMatrix

#pragma once
#include "AbstractRealMatrix.hpp"

class BandMatrix : public AbstractRealMatrix { // https://en.wikipedia.org/wiki/Band_matrix
	size_t _w, _p; // w := width of band (i.e. numb of non-zero diagonals), p := (w – 1) / 2
	bool _isDecomposed; // this flag shows whether matrix is LU decomposed 
	double* _A, ** _L, * _D, ** _U; // A = [ L | D | U ] := n × w band matrix (lower triangle, diagonal, upper triangle)
	bool _computeLU();
	void _forwardSubst(std::vector<double>&) const; // forward substitution; find vector y such that L.y = b
	void _backSubst(std::vector<double>&) const; // back substitution; find vector x such that U.x = y
public:
	BandMatrix(size_t, size_t); // create (allocate memory for) n × w matrix
	~BandMatrix(); // free memory
	// virtual methods to be implemented
	double& set(size_t, size_t);
	double get(size_t, size_t) const;
	std::vector<double> solve(std::vector<double> const &);
	std::istream& load(std::istream&);
	std::ostream& save(std::ostream&) const;
	std::vector<double> mult(std::vector<double> const &) const;
};