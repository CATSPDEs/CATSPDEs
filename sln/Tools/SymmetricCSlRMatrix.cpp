#include <complex>
#include "SymmetricCSlRMatrix.hpp"

template <typename T>
SymmetricCSlRMatrix::SymmetricCSlRMatrix(size_t n, size_t nnz)
	: AbstractSparseMatrix(n, n)
	, _diag(n)
	, _iptr(n + 1)
	, _jptr(nnz)
	, _lval(nnz) {
	_iptr[n] = nnz;
}

template <typename T>
SymmetricCSlRMatrix::SymmetricCSlRMatrix(AdjacencyList const & adjList) 
	: AbstractSparseMatrix(adjList.size(), adjList.size())
	, _iptr(adjList.size() + 1)
	, _diag(adjList.size()) {
	for (size_t i = 0; i < adjList.size(); i++)
		_iptr[i + 1] = _iptr[i] + adjList[i].size();
	_jptr.reserve(_iptr[_w]);
	for (size_t i = 0; i < adjList.size(); i++)
		for (auto neighbour : adjList[i])
			_jptr.emplace_back(neighbour);
	_lval.resize(_iptr[_w]);
}

// private methods

template <typename T>
T& SymmetricCSlRMatrix::_set(size_t i, size_t j) {
	if (i == j) return _diag[i];
	if (i < j) swap(i, j); // lower triangular part
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _lval[k];
		else if (_jptr[k] > j) throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

template <typename T>
T SymmetricCSlRMatrix::_get(size_t i, size_t j) const {
	if (i == j) return _diag[i];
	if (i < j) swap(i, j);
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _lval[k];
		else if (_jptr[k] > j) return 0.;
	return 0.;
}

// public methods

template <typename T>
istream & SymmetricCSlRMatrix::loadSparse(istream& input) {
	// we know matrix size (_w) and numb of nonzeroes 
	// in lower triangular part (_iptr[_w])
	// so here is istream structure:
	//
	// (1) _w + 1    elements of _iptr (the last element is already known),
	// (2) _iptr[_w] elements of _jptr,
	// (3) _iptr[_w] elements of _lptr,
	// (4) _w        elements of _diag
	//
	size_t i, max = _w + 1;
	for (i = 0; i < max; ++i)
		input >> _iptr[i];
	max = _iptr[_w];
	for (i = 0; i < max; ++i)
		input >> _jptr[i];
	for (i = 0; i < max; ++i)
		input >> _lval[i];
	max = _w;
	for (i = 0; i < max; ++i)
		input >> _diag[i];
	return input;
}

template <typename T>
ostream& SymmetricCSlRMatrix::saveSparse(ostream& output) const {
	return output << _w << ' ' << _iptr[_w] << '\n'
		<< _iptr << '\n' 
		<< _jptr << '\n'
		<< _lval << '\n'
		<< _diag;
}

template <typename T>
SymmetricCSlRMatrix SymmetricCSlRMatrix::ILDL()
{
	SymmetricCSlRMatrix L(*this);
	for (size_t k = 1;k < _w;k++){
		for (size_t j = _iptr[k]; j < _iptr[k + 1];j++)	{
			for (size_t i = _iptr[k];i < j;i++)
				L._lval[j] -= L._lval[i] * L._diag[_jptr[i]] *
				[&](int q, int d) {
				int p;	
				for (p = _iptr[d]; p < _iptr[d + 1];p++)
						if (_jptr[p] == q) return L._lval[p];
					return 0.;
				}(_jptr[i], _jptr[j]);
			L._lval[j]/= L._diag[_jptr[j]];
		}
		for (size_t j = _iptr[k];j < _iptr[k + 1];j++)
			L._diag[k] -= L._lval[j] * L._lval[j] * L._diag[_jptr[j]];
	}
	return L;
}

template <typename T>
vector<T> SymmetricCSlRMatrix::_mult(vector<T> const &u) const { // compute v := A.u
	vector<T> v(_w, 0.);
	for (size_t i = 0; i < _w; ++i) {
		v[i] += u[i] * _diag[i];
		for (size_t j = _iptr[i]; j < _iptr[i + 1]; ++j) {
			v[i]        += _lval[j] * u[_jptr[j]];
			v[_jptr[j]] += _lval[j] * u[i];
		}
	}
	return v;
}

template <typename T>
SymmetricCSlRMatrix& SymmetricCSlRMatrix::setZero() {
	fill(_diag.begin(), _diag.end(), 0.);
	fill(_lval.begin(), _lval.end(), 0.);
	return *this;
}


template <typename T>
vector<T> SymmetricCSlRMatrix::forwardSubstitution(vector<T> const &f) const
{
	vector<T> res(_w);
	vector<T> temp(f);
	for (int i = 0;i < _w;i++) {
		for (int j = _iptr[i];j < _iptr[i + 1];j++)
			temp[i] = f[i] - res[_jptr[j]] * _lval[j];
		res[i] = temp[i] / _diag[i];
	}
	return res;
}

template <typename T>
vector<T> SymmetricCSlRMatrix::backwardSubstitution(vector<T> const &f) const
{
	return vector<T>(_w, 0);
}

template <typename T>
SymmetricCSlRMatrix& SymmetricCSlRMatrix::decompose() {
	// ... 
	return *this;
}

// explicit instantiation

template class SymmetricCSlRMatrix<double>;
template class SymmetricCSlRMatrix<complex<double>>;
