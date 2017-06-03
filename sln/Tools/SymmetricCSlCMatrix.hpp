﻿#pragma once
#include <type_traits> // is_same 
#include "CSlCMatrix.hpp"

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
	, public AbstractPreconditioner<T> {
	// fancy name, yeah
	// stands for Symmetric Compressed Sparse (lower triangular) Column
	std::vector<Index> _colptr, // vector of column pointers and
	                   _rowind; // ″         row indicies (see example to get what these guys are)
	std::vector<T>     _lval, // vector of elements of lower triangular part of matrix (raw by raw)
	                   _diag; // ″         diagonal elements (for FEM / FVM is always > 0, 
			                   // so we store it explicitly)
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
	// _rowind = { 2, 4, 2, 4, 3, 5, 4, 5 }
	// _colptr = { 0, 2, 4, 6, 7, 8, 8 }
	// 
	// so, for one, 1st col (we count from zero) countains _colptr[1] – _colptr[0] = 2 – 0 = 2 elements, 
	// starting w/ element _lval[_colptr[1]] = _lval[2] = 5 in (_rowind[_colptr[1]] = _rowind[2] = 2)nd row
	// size of matrix is 6, so _colptr[6] = 8 is numb of elements in _lval or _rowind
	// makes sense!
	//
	// virtual methods to be implemented
	T& _set(Index, Index) final;
	T  _get(Index, Index) const final;
public:
	SymmetricCSlCMatrix(std::vector<Index> const & colptr, std::vector<Index> const & rowind, std::vector<T> const & lval, std::vector<T> const & diag) : AbstractMatrix(diag.size()), _colptr(colptr), _rowind(rowind), _lval(lval), _diag(diag) {}
	explicit SymmetricCSlCMatrix(Index n = 1, Index nnz = 0); // order of matrix and numb of nonzero elems in lower triangular part
	// virtual methods to be implemented
	Index nnz() const final { return _colptr[_w] + _w; } // “+ _w” because of _diag
	SymmetricCSlCMatrix& operator=(T const &) final;
	void mult(T const * by, T* result) const final;
	// i/o
	SymmetricCSlCMatrix& importSparse(std::istream& from = cin) final;
	void                 exportSparse(std::ostream& to = cout) const final;
	// almost same methods from base class (in order to work w/ strings instead of streams)
	using AbstractSparseMatrix::importSparse;
	using AbstractSparseMatrix::exportSparse;
	HarwellBoeingHeader importHarwellBoeing(std::string const &) final; // load from Harwell–Boeing file
	void                exportHarwellBoeing(std::string const &, Parameters const & params = {}) const final;
	SymmetricCSlCMatrix& generatePatternFrom(DOFsConnectivityList const &) final;
	SymmetricCSlCMatrix& enforceDirichletBCs(Index2Value<T> const &, T*) final;
	// precond
	std::vector<T> forwSubst(std::vector<T> const &, double w = 1.) const final;
	std::vector<T> backSubst(std::vector<T> const &, double w = 1.) const final;
	std::vector<T> diagSubst(std::vector<T> const &) const final;
	std::vector<T> multDiag(std::vector<T> const &) const final;
	SymmetricCSlCMatrix<T>& decompose() final;
	// explicit conversion to non-symmetric format http://en.cppreference.com/w/cpp/language/cast_operator
	explicit operator CSlCMatrix<T>() const {
		return { _colptr, _rowind, _lval, _lval, _diag };
	}
};

// implementation

template <typename T>
SymmetricCSlCMatrix<T>::SymmetricCSlCMatrix(Index n, Index nnz)
	: AbstractMatrix<T>(n)
	, _colptr(n + 1)
	, _rowind(nnz)
	, _lval(nnz)
	, _diag(n) {
	_colptr[n] = nnz;
}

// private methods

template <typename T>
T& SymmetricCSlCMatrix<T>::_set(Index i, Index j) {
	if (i == j) return _diag[i];
	if (i < j) std::swap(i, j); // switch to lower triangular part
	for (Index k = _colptr[j]; k < _colptr[j + 1]; ++k)
		if (_rowind[k] == i) return _lval[k];
	throw std::invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
}

template <typename T>
T SymmetricCSlCMatrix<T>::_get(Index i, Index j) const {
	if (i == j) return _diag[i];
	if (i < j) std::swap(i, j);
	for (Index k = _colptr[j]; k < _colptr[j + 1]; ++k)
		if (_rowind[k] == i) return _lval[k];
	return 0.;
}

// public methods

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::operator=(T const & val) {
	std::fill(_lval.begin(), _lval.end(), val);
	std::fill(_diag.begin(), _diag.end(), val);
	return *this;
}

template <typename T>
void SymmetricCSlCMatrix<T>::mult(T const * by, T* result) const { // = multByTranspose() since A = A^T
	for (Index j = 0; j < _w; ++j) {
		result[j] += _diag[j] * by[j];
		for (Index i = _colptr[j]; i < _colptr[j + 1]; ++i) {
			result[_rowind[i]] += _lval[i] * by[j];
			result[j]          += _lval[i] * by[_rowind[i]];
		}
	}
}

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::importSparse(std::istream& input) {
	// stdin structure:
	// (1) n =: _h = _w (size of square matrix)
	// (2) numb of nonzeros in lower triangular part =: _colptr[n + 1]
	// (3) n + 1    elements of _colptr,
	// (4) _colptr[n] elements of _rowind,
	// (5) _colptr[n] elements of _lval, and
	// (6) n        elements of _diag
	Index nnz;
	input >> _h >> nnz;
	_w = _h;
	_colptr.resize(_w + 1);
	_rowind.resize(nnz);
	_lval.resize(nnz);
	_diag.resize(_w);
	input >> _colptr >> _rowind >> _lval >> _diag;
	return *this;
}

template <typename T>
void SymmetricCSlCMatrix<T>::exportSparse(std::ostream& output) const {
	output << _w << ' ' << _colptr[_w] << '\n'
	       << _colptr << '\n'
	       << _rowind << '\n'
	       << _lval << '\n'
	       << _diag;
}

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::generatePatternFrom(DOFsConnectivityList const & list) {
	// (1) compute column pointers
	for (Index i = 0; i < _w; ++i) _colptr[i + 1] = _colptr[i] + list[i].size();
	// (2) compute row indicies
	_rowind.reserve(_colptr[_w]);
	_lval.resize(_colptr[_w]);
	for (Index col = 0; col < _w; ++col)
		for (auto row : list[col])
			_rowind.emplace_back(row);
	return *this;
}

template <typename T>
HarwellBoeingHeader SymmetricCSlCMatrix<T>::importHarwellBoeing(std::string const & fileName) {
	throw std::logic_error("not implemented");
	// load header
	HarwellBoeingHeader header;
	loadHarwellBoeingHeader_f90(fileName.c_str(), &header);
	// check matrix type
	if (header.mxtype[1] != 'S') throw std::invalid_argument("SymmetricCSlC type is useful for symmetric (.*s*) Harwell-Boeing matrices; try CSC for (.*r*) or (.*u*) types");
	if (header.mxtype[2] != 'A') throw std::invalid_argument("only assembled (.**a) Harwell-Boeing matrices are supported");
	if (std::is_same<T, double>::value) { // real matrix 
		if (header.mxtype[0] == 'C') throw std::invalid_argument("instantiate matrix w/ complex<double> in order to work w/ complex (.c**) Harwell-Boeing matrices");
	}
	else // complex matrix
		if (header.mxtype[0] == 'R') throw std::invalid_argument("instantiate matrix w/ double in order to work w/ real (.r**) Harwell-Boeing matrices");
	
	// read structure and values
	//// HB stores diagonal in lower triangular part, but we in CATSPDEs prefer to store it
	//// separately since in FEM / FDM / FVM diags are nonzero
	//// so we need to fix column pointers, row indicies, and separate diagonal from values array
	//vector<Index> rowind(_colptr[_w]);
	//vector<T>      values(_colptr[_w]);
	//Index diagSize; // numb of diag elements stored as a part of lower triangle of HB matrix
	//loadHarwellBoeingStruct_f90(
	//	fileName.c_str(), 
	//	&header,
	//	_colptr.data(), 
	//	rowind.data(),
	//	reinterpret_cast<double*>(values.data()), // http://en.cppreference.com/w/cpp/numeric/complex > non-static data members
	//);

	//for (i = )

	//// shrink
	//_rowind.resize(_rowind.size() - diagSize);
	//_lval.resize(_lval.size() - diagSize);
	//for_each(_colptr.begin(), _colptr.end(), [](Index& i) { --i; });
	//for_each(_rowind.begin(), _rowind.end(), [](Index& i) { --i; });

	return header;
}

template <typename T>
void SymmetricCSlCMatrix<T>::exportHarwellBoeing(std::string const & fileName, Parameters const & params) const {
	throw std::logic_error("not implemented");
}

template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::enforceDirichletBCs(Index2Value<T> const & ind2val, T* rhs) {
	for (Index i = 0; i < getOrder(); ++i) { // for each row (column)…
		auto kvpIter = ind2val.find(i);
		if (kvpIter != ind2val.end()) {
			for (Index j = _colptr[i]; j < _colptr[i + 1]; ++j) {
				rhs[_rowind[j]] -= _lval[j] * kvpIter->second;
				_lval[j] = 0.;
			}
			_diag[i] = 1.;
			rhs[i] = kvpIter->second;
		}
		else
			for (Index j = _colptr[i]; j < _colptr[i + 1]; ++j)
				if ((kvpIter = ind2val.find(_rowind[j])) != ind2val.end()) {
					rhs[i] -= _lval[j] * kvpIter->second;
					_lval[j] = 0.;
				}
	}
	return *this;
}

// find y := [ L + wD ]^-1 . x
template <typename T>
std::vector<T> SymmetricCSlCMatrix<T>::forwSubst(std::vector<T> const & x, double w = 1.) const {
	std::vector<T> y { x };
	Index i, j, k;
	for (j = 0; j < getOrder(); ++j) {
		if (w) y[j] /= (w * _diag[j]);
		for (k = _colptr[j]; k < _colptr[j + 1]; ++k) {
			i = _rowind[k];
			auto& l_ij = _lval[k];
			y[i] -= l_ij * y[j];
		}
	}
	return y;
}

// find y := [ U + wD ]^-1 . x
template <typename T>
std::vector<T> SymmetricCSlCMatrix<T>::backSubst(std::vector<T> const & x, double w = 1.) const {
	std::vector<T> y { x };
	SignedIndex i;
	Index j, k;
	for (i = getOrder() - 1; i >= 0; --i) {
		for (k = _colptr[i]; k < _colptr[i + 1]; ++k) {
			j = _rowind[k];
			auto& u_ij = _lval[k]; // = l_ji
			y[i] -= u_ij * y[j];
		}
		if (w) y[i] /= (w * _diag[i]);
	}
	return y;
}

template <typename T>
std::vector<T> SymmetricCSlCMatrix<T>::diagSubst(std::vector<T> const & x) const {
	std::vector<T> y(getOrder());
	for (Index i = 0; i < getOrder(); ++i)
		y[i] = x[i] / _diag[i];
	return y;
}

template <typename T>
std::vector<T> SymmetricCSlCMatrix<T>::multDiag(std::vector<T> const & x) const {
	std::vector<T> y(getOrder());
	for (Index i = 0; i < getOrder(); ++i)
		y[i] = x[i] * _diag[i];
	return y;
}

// ILDL^T(0) Crout decomposition from Saad, p. 333
// A ~ (L + I) . D . (I + L^T) 
template <typename T>
SymmetricCSlCMatrix<T>& SymmetricCSlCMatrix<T>::decompose() {
	Index n = getOrder();
	auto colptr_first = _colptr;
	std::vector<std::list<Index>> nnz_rows_of_col(n);
	std::vector<T> kth_col(n);
	for (Index k = 0; k < n; ++k) {
		// copy column and row to be updated on kth step
		for (Index m = _colptr[k]; m < _colptr[k + 1]; ++m) {
			Index i = _rowind[m];
			kth_col[i] = _lval[m];
		}
		// loop over only those cols which get multiplied by a nnz row
		for (Index j : nnz_rows_of_col[k]) {
			auto l_kj = _get(k, j);
			_diag[k] -= _diag[j] * l_kj * l_kj;
			while (colptr_first[j] < _colptr[j + 1] && _rowind[colptr_first[j]] <= k) ++colptr_first[j];
			for (Index m = colptr_first[j]; m < _colptr[j + 1]; ++m) {
				Index i = _rowind[m]; auto l_ij = _lval[m];
				kth_col[i] -= _diag[j] * l_kj * l_ij;
			}
		}
		// update nnz_rows_of_col lists, kth column, and kth row of the matrix
		for (Index m = _colptr[k]; m < _colptr[k + 1]; ++m) {
			Index i = _rowind[m];
			nnz_rows_of_col[i].push_back(k);
			_lval[m] = kth_col[i] / _diag[k];
		}
	}
	return *this;
}