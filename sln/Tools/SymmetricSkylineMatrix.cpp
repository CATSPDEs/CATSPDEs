#include "SymmetricSkylineMatrix.hpp"

std::istream & SymmetricSkylineMatrix::loadSparse(std::istream& input) {
	// we know matrix size (_n) and numb of nonzeroes 
	// in lower triangular part (_iptr[_n])
	// so here is istream structure:
	//
	// (1) _n - 1    elements of _iptr, w/o _iptr[_n]
	// (2) _iptr[_n] elements of _jptr,
	// (3) _iptr[_n] elements of _lptr,
	// (4) _n        elements of _diag
	//
	size_t i, max = _n;
	for (i = 0; i < max; ++i)
		input >> _iptr[i];
	max = _iptr[_n];
	for (i = 0; i < max; ++i)
		input >> _jptr[i];
	for (i = 0; i < max; ++i)
		input >> _lval[i];
	max = _n;
	for (i = 0; i < max; ++i)
		input >> _diag[i];
	return input;
}

std::ostream& SymmetricSkylineMatrix::saveSparse(std::ostream& output) const
{
	return output << _n << ' ' << _iptr[_n] << '\n'
		<< std::vector<size_t>(_iptr.begin(), _iptr.end() - 1) << '\n' 
		// we do not want last element to be written—we already have it
		<< _jptr << '\n'
		<< _lval << '\n'
		<< _diag;
}

double & SymmetricSkylineMatrix::_set(size_t, size_t)
{
	//placeholder
	return _lval[0];
}

double SymmetricSkylineMatrix::_get(size_t, size_t) const
{
	//placeholder
	return 0.0;
}

SymmetricSkylineMatrix::SymmetricSkylineMatrix(size_t n, size_t nnz)
	: AbstractSparseMatrix(n)
	, _iptr(n + 1)
	, _jptr(nnz)
	, _lval(nnz) {
	_iptr[n] = nnz;
}

SymmetricSkylineMatrix::SymmetricSkylineMatrix(AdjacencyList adjList):
	AbstractSparseMatrix(adjList.size()),
	_iptr(adjList.size()+1),
	_diag(adjList.size()){
	size_t nnz = 0;
	_jptr.reserve((_n*_n)/2);
	for (size_t i = 0;i < adjList.size();i++) {
		_iptr[i + 1] = _iptr[i] + adjList[i].size();
		for (auto neighbour : adjList[i]) 
			_jptr.push_back(neighbour);
	}
	_lval.resize(_iptr[_n]);
	_jptr.shrink_to_fit();
}

std::vector<double> SymmetricSkylineMatrix::solve(std::vector<double> const &)
{
	//placeholder
	return std::vector<double>();
}

std::vector<double> SymmetricSkylineMatrix::mult(std::vector<double> const &) const
{
	//placeholder
	return std::vector<double>();
}
