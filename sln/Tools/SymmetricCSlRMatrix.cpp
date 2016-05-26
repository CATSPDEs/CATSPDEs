#include "SymmetricCSlRMatrix.hpp"

SymmetricCSlRMatrix::SymmetricCSlRMatrix(size_t n, size_t nnz)
	: AbstractSparseMatrix(n)
	, _iptr(n + 1)
	, _jptr(nnz)
	, _lval(nnz) {
	_iptr[n] = nnz;
}

SymmetricCSlRMatrix::SymmetricCSlRMatrix(AdjacencyList const & adjList) 
	: AbstractSparseMatrix(adjList.size())
	, _iptr(adjList.size() + 1)
	, _diag(adjList.size()) {
	for (size_t i = 0;i < adjList.size();i++)
		_iptr[i + 1] = _iptr[i] + adjList[i].size();
	_jptr.reserve(_iptr[_n]);
	for (size_t i = 0;i < adjList.size();i++)
		for (auto neighbour : adjList[i])
			_jptr.emplace_back(neighbour);
	_lval.resize(_iptr[_n]);
}

// private methods

double& SymmetricCSlRMatrix::_set(size_t i, size_t j) {
	if (i == j) return _diag[i];
	if (i < j) std::swap(i, j); // lower triangular part
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _lval[k];
		else if (_jptr[k] > j) throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

double SymmetricCSlRMatrix::_get(size_t i, size_t j) const {
	if (i == j) return _diag[i];
	if (i < j) std::swap(i, j);
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _lval[k];
		else if (_jptr[k] > j) return 0.;
}

// public methods

std::istream & SymmetricCSlRMatrix::loadSparse(std::istream& input) {
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

std::ostream& SymmetricCSlRMatrix::saveSparse(std::ostream& output) const {
	return output << _n << ' ' << _iptr[_n] << '\n'
		<< std::vector<size_t>(_iptr.begin(), _iptr.end() - 1) << '\n' 
		// we do not want last element to be written—we already have it
		<< _jptr << '\n'
		<< _lval << '\n'
		<< _diag;
}

std::vector<double> SymmetricCSlRMatrix::solve(std::vector<double> const &) {
	//placeholder
	return std::vector<double>();
}

std::vector<double> SymmetricCSlRMatrix::mult(std::vector<double> const &u) const { // compute v := A.u
	std::vector<double> v(_n, 0.);
	for (size_t i = 0; i < _n; ++i) {
		v[i] += u[i] * _diag[i];
		for (size_t j = _iptr[i]; j < _iptr[i + 1]; ++j) {
			v[i]        += _lval[j] * u[_jptr[j]];
			v[_jptr[j]] += _lval[j] * u[i];
		}
	}
	return v;
}
