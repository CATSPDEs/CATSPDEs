/*
	Alexander Žilyakov, Sep 2016
*/
#pragma once
#include <type_traits> // is_same 
#include "AbstractHarwellBoeingMatrix.hpp"
#include "AbstractMultipliableMatrix.hpp"
// TODO:
//#include "IDecomposable.hpp"
//#include "AdjacencyList.hpp"

template <typename T> class SymmetricCSlCMatrix;
// matrix is stored in symmetric CSlC format…

template <typename T> using SymmetricHBMatrix = SymmetricCSlCMatrix<T>;
// …aka Harwell–Boeing format…

template <typename T>
class SymmetricCSlCMatrix
	: public AbstractHarwellBoeingMatrix<T>
	, public AbstractMultipliableMatrix<T> {
	// fancy name, yeah
	// stands for Symmetric Compressed Sparse (lower triangular) Column
	std::vector<double> _lval, // vector of elements of lower triangular part of matrix (raw by raw)
	                    _diag; // ″         diagonal elements (for FEM / FVM is always > 0, 
	                           //                              so we store it explicitly)
	std::vector<size_t> _jptr, // vector of column pointers and
	                    _iptr; // ″         row indicies (see example to get what these guys are)
	// example of 6 × 6 matrix:
	//		
	//    8 ─ ─ ─ ─ ┐
	//    × 1   sym |
	//    2 5 4     |
	//    × × 3 2   |
	//    7 8 × 1 7 |
	//    × × 2 × 4 2
	//
	// (“×” denotes zero element we do not store)
	//
	// _diag = { 8, 1, 4, 2, 7, 2 }
	// _lval = { 2, 7, 5, 8, 3, 2, 1, 4 }
	// _iptr = { 2, 4, 2, 4, 3, 5, 4, 5 }
	// _jptr = { 0, 2, 4, 6, 7, 8, 8 }
	// 
	// so, for one, 1st col (we count from zero) countains _jptr[1] – _jptr[0] = 2 – 0 = 2 elements, 
	// starting w/ element _lval[_jptr[1]] = _lval[2] = 5 in (_iptr[_jptr[1]] = _iptr[2] = 2)nd row
	// size of matrix is 6, so _jptr[6] = 8 is numb of elements in _lval or _iptr
	// makes sense!
	//
	// virtual methods to be implemented
	size_t _nnz() const final { return _jptr[_w] + _w; } // “+ _w” because of _diag
	T& _set(size_t, size_t) final;
	T  _get(size_t, size_t) const final;
public:
	SymmetricCSlCMatrix(size_t n, size_t nnz); // order of matrix and numb of nonzero elems in lower triangular part
	// TODO:
	//SymmetricCSlCMatrix(AdjacencyList const &); // generate matrix portrait from adjacency list of mesh nodes 
	~SymmetricCSlRMatrix() {}
	// virtual methods to be implemented
	SymmetricCSlCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	SymmetricCSlCMatrix& loadSparse(istream& from = cin) final;
	void                 saveSparse(ostream& to = cout) const final;
	SymmetricCSlCMatrix& loadHarwellBoeing(HarwellBoeingHeader const &, string const &) final; // load from Harwell–Boeing file
	void                 saveHarwellBoeing(string const &, Parameters const & params = {}) const final;
	// TODO:
	//SymmetricCSlCMatrix& decompose() override;
	//vector<T> forwardSubstitution(vector<T> const &) const override;
	//vector<T> backwardSubstitution(vector<T> const &) const override;
};

// implementation

template <typename T>
SymmetricCSlCMatrix<T>::SymmetricCSlCMatrix(size_t n, size_t nnz)
	: AbstractMatrix<T>(n, n)
	, _jptr(n + 1)
	, _iptr(nnz)
	, _lval(nnz)
	, _diag(n) {
	_jptr[n] = nnz;
}

// private methods

template <typename T>
T& SymmetricCSlCMatrix<T>::_set(size_t i, size_t j) {
	if (i == j) return _diag[i];
	if (i < j) swap(i, j); // switch to lower triangular part
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _lval[k];
		else if (_iptr[k] > i) throw invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
	throw invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
}

template <typename T>
T SymmetricCSlCMatrix<T>::_get(size_t i, size_t j) const {
	if (i == j) return _diag[i];
	if (i < j) swap(i, j);
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _lval[k];
		else if (_iptr[k] > i) return 0.;
	return 0.;
}

// public methods

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::operator=(T const & val) {
	fill(_mval.begin(), _lval.end(), val);
	fill(_diag.begin(), _diag.end(), val);
	return *this;
}

template <typename T>
void SymmetricCSlCMatrix<T>::mult(T const * by, T* result) const { // = multByTranspose() since A = A^T
	size_t i, j;
	for (j = 0; j < _w; ++j) {
		result[j] = _diag[j] * by[j];
		for (i = _jptr[j]; i < _jptr[j + 1]; ++i)
			result[j] += _lval[i] * by[_iptr[i]];
	}
}

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::loadSparse(istream& input) {
	// stdin structure:
	// (1) _w + 1    elements of _jptr,
	// (2) _jptr[_w] elements of _iptr,
	// (3) _jptr[_w] elements of _lval, and
	// (4) _w        elements of _diag
	input >> _jptr >> _iptr >> _lval >> _diag;
	return *this;
}

template <typename T>
void SymmetricCSlCMatrix<T>::saveSparse(ostream& output) const {
	output << _w << ' ' << _jptr[_w] << '\n'
	       << _jptr << '\n'
	       << _iptr << '\n'
	       << _lval << '\n'
	       << _diag;
}

extern "C" void loadHarwellBoeingStructSym_f90(char const *, HarwellBoeingHeader const *, size_t*, size_t*, double*); 
// read struct and values—almost like loadHarwellBoeingStruct_f90(…) routine—but for symmetric matrices
// HB stores diagonal in lower triangular part, but we in CATSPDEs prefer to store it
// separately since in FEM / FDM / FVM diags are nonzero
// so we need slightly modified version of the loadHarwellBoeingStruct_f90(…) func for symmetric matrices
template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::loadHarwellBoeing(HarwellBoeingHeader const & header, string const & fileName) {
	// check matrix type
	if (header.mxtype[1] != 'S') throw invalid_argument("SymmetricCSlC type is useful for symmetric (.*s*) Harwell-Boeing matrices; try CSC for (.*r*) or (.*u*) types");
	if (header.mxtype[2] != 'A') throw invalid_argument("only assembled (.**a) Harwell-Boeing matrices are supported");
	if (is_same<T, double>::value) { // real matrix 
		if (header.mxtype[0] == 'C') throw invalid_argument("instantiate matrix w/ complex<double> in order to work w/ complex (.c**) Harwell-Boeing matrices");
	}
	else // complex matrix
		if (header.mxtype[0] == 'R') throw invalid_argument("instantiate matrix w/ double in order to work w/ real (.r**) Harwell-Boeing matrices");
	// read structure and values
	loadHarwellBoeingStructSym_f90(fileName.c_str(), &header, _jptr.data(), _iptr.data(),
	                               reinterpret_cast<double*>(_lval.data())); // http://en.cppreference.com/w/cpp/numeric/complex > non-static data members
	// TODO: FIX _iptr, _lval sizes!!
	for_each(_jptr.begin(), _jptr.end(), [](size_t& i) { --i; });
	for_each(_iptr.begin(), _iptr.end(), [](size_t& i) { --i; });
	return *this;
}

extern "C" void saveHarwellBoeingStructSym_f90(char const *, HarwellBoeingHeader const *, size_t const *, size_t const *, double const *);
template <typename T>
void SymmetricCSlCMatrix<T>::saveHarwellBoeing(string const & fileName, Parameters const & params) const {
	// …
}