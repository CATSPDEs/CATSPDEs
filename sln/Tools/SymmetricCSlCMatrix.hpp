#pragma once
#include <type_traits> // is_same 
#include "AbstractFEMatrix.hpp"
#include "AbstractHarwellBoeingMatrix.hpp"
#include "AbstractMultipliableMatrix.hpp"
// TODO:
//#include "IDecomposable.hpp"

/*
	Alexander Žilyakov, Sep 2016
*/

template <typename T> class SymmetricCSlCMatrix;
// matrix is stored in symmetric CSlC format…

template <typename T> using SymmetricHBMatrix = SymmetricCSlCMatrix<T>;
// …aka Harwell–Boeing format…

template <typename T>
class SymmetricCSlCMatrix
	: public AbstractFEMatrix<T>
	, public AbstractHarwellBoeingMatrix<T>
	, public AbstractMultipliableMatrix<T> {
	// fancy name, yeah
	// stands for Symmetric Compressed Sparse (lower triangular) Column
	vector<size_t> _jptr, // vector of column pointers and
	               _iptr; // ″         row indicies (see example to get what these guys are)
	vector<double> _lval, // vector of elements of lower triangular part of matrix (raw by raw)
	               _diag; // ″         diagonal elements (for FEM / FVM is always > 0, 
	                      //                              so we store it explicitly)
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
	explicit SymmetricCSlCMatrix(size_t n = 1, size_t nnz = 0); // order of matrix and numb of nonzero elems in lower triangular part
	// TODO: 
	~SymmetricCSlCMatrix() {}
	// virtual methods to be implemented
	SymmetricCSlCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	SymmetricCSlCMatrix& loadSparse(istream& from = cin) final;
	void                 saveSparse(ostream& to = cout) const final;
	SymmetricCSlCMatrix& generatePatternFrom(AdjacencyList const &) final;
	SymmetricCSlCMatrix& loadHarwellBoeing(string const &, HarwellBoeingHeader* headerPtr = nullptr) final; // load from Harwell–Boeing file
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
	fill(_lval.begin(), _lval.end(), val);
	fill(_diag.begin(), _diag.end(), val);
	return *this;
}

template <typename T>
void SymmetricCSlCMatrix<T>::mult(T const * by, T* result) const { // = multByTranspose() since A = A^T
	size_t i, j;
	for (j = 0; j < _w; ++j) result[j] = _diag[j] * by[j];
	for (j = 0; j < _w; ++j) 
		for (i = _jptr[j]; i < _jptr[j + 1]; ++i) {
			result[_iptr[i]] += _lval[i] * by[j];
			result[j]        += _lval[i] * by[_iptr[i]];
		}
}

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::loadSparse(istream& input) {
	// stdin structure:
	// (1) n =: _h = _w (size of square matrix)
	// (2) numb of nonzeros in lower triangular part =: _jptr[n + 1]
	// (3) n + 1    elements of _jptr,
	// (4) _jptr[n] elements of _iptr,
	// (5) _jptr[n] elements of _lval, and
	// (6) n        elements of _diag
	size_t nnz;
	input >> _h >> nnz;
	_w = _h;
	_jptr.resize(_w + 1);
	_iptr.resize(nnz);
	_lval.resize(nnz);
	_diag.resize(_w);
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

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::generatePatternFrom(AdjacencyList const & adjList) {
	// TODO: check!
	_h = _w = adjList.size();
	_jptr.resize(_w + 1);
	_diag.resize(_w);
	for (size_t j = 0; j < _w; ++j)
		_jptr[j + 1] = _jptr[j] + adjList[j].size();
	_iptr.reserve(_jptr[_w]);
	for (size_t i = 0; i < adjList.size(); i++)
		for (auto neighbour : adjList[i])
			_iptr.emplace_back(neighbour);
	_lval.resize(_jptr[_w]);
	return *this;
}

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::loadHarwellBoeing(string const & fileName, HarwellBoeingHeader* headerPtr) {
	//// check matrix type
	//if (header.mxtype[1] != 'S') throw invalid_argument("SymmetricCSlC type is useful for symmetric (.*s*) Harwell-Boeing matrices; try CSC for (.*r*) or (.*u*) types");
	//if (header.mxtype[2] != 'A') throw invalid_argument("only assembled (.**a) Harwell-Boeing matrices are supported");
	//if (is_same<T, double>::value) { // real matrix 
	//	if (header.mxtype[0] == 'C') throw invalid_argument("instantiate matrix w/ complex<double> in order to work w/ complex (.c**) Harwell-Boeing matrices");
	//}
	//else // complex matrix
	//	if (header.mxtype[0] == 'R') throw invalid_argument("instantiate matrix w/ double in order to work w/ real (.r**) Harwell-Boeing matrices");
	//// read structure and values
	//// HB stores diagonal in lower triangular part, but we in CATSPDEs prefer to store it
	//// separately since in FEM / FDM / FVM diags are nonzero
	//// so we need to fix column pointers, row indicies, and separate diagonal from values array
	//vector<size_t> rowind(_jptr[_w]);
	//vector<T>      values(_jptr[_w]);
	//size_t diagSize; // numb of diag elements stored as a part of lower triangle of HB matrix
	//loadHarwellBoeingStruct_f90(
	//	fileName.c_str(), 
	//	&header,
	//	_jptr.data(), 
	//	rowind.data(),
	//	reinterpret_cast<double*>(values.data()), // http://en.cppreference.com/w/cpp/numeric/complex > non-static data members
	//);

	//for (i = )

	//// shrink
	//_iptr.resize(_iptr.size() - diagSize);
	//_lval.resize(_lval.size() - diagSize);
	//for_each(_jptr.begin(), _jptr.end(), [](size_t& i) { --i; });
	//for_each(_iptr.begin(), _iptr.end(), [](size_t& i) { --i; });
	return *this;
}

template <typename T>
void SymmetricCSlCMatrix<T>::saveHarwellBoeing(string const & fileName, Parameters const & params) const {
	// …
}