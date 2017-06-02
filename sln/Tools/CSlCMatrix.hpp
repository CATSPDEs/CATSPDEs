#pragma once
#include "AbstractPreconditioner.hpp"
#include "CSCMatrix.hpp"

/*
	Alexander Žilyakov, Jan 2017
*/

template <typename T>
class CSlCMatrix
	: public AbstractFEMatrix<T>
	, public AbstractPreconditioner<T> {
	std::vector<Index> _colptr, _rowind; 
	std::vector<T>     _lval, _uval, _diag; 
	T& _set(Index, Index) final;
	T  _get(Index, Index) const final;
public:
	CSlCMatrix(std::vector<Index> const & colptr, std::vector<Index> const & rowind, std::vector<T> const & lval, std::vector<T> const & uval, std::vector<T> const & diag) : AbstractMatrix(diag.size()), _colptr(colptr), _rowind(rowind), _lval(lval), _uval(uval), _diag(diag) {}
	explicit CSlCMatrix(Index n = 1, Index nnz = 0); // order of matrix and numb of nonzero elems in lower (upper) triangular part
	// virtual methods to be implemented
	Index nnz() const final { return 2 * _colptr[_w] + _w; } // “+ _w” because of _diag
	CSlCMatrix& operator=(T const &) final;
	void mult(T const * by, T* result) const final;
	// i/o
	CSlCMatrix& importSparse(std::istream& from = cin) final;
	void        exportSparse(std::ostream& to = cout) const final;
	// almost same methods from base class (in order to work w/ strings instead of streams)
	using AbstractSparseMatrix::importSparse;
	using AbstractSparseMatrix::exportSparse;
	CSlCMatrix& generatePatternFrom(DOFsConnectivityList const &) final;
	CSlCMatrix& enforceDirichletBCs(Index2Value<T> const &, T*) final;
	// precond
	std::vector<T> forwSubst(std::vector<T> const & x, double w = 1.) const final;
	std::vector<T> backSubst(std::vector<T> const &, double w = 1.) const final;
	std::vector<T> diagSubst(std::vector<T> const & x) const final;
	std::vector<T> multDiag(std::vector<T> const & x) const final;
	CSlCMatrix<T>& decompose() final;
	// explicit conversion to non-symmetric format http://en.cppreference.com/w/cpp/language/cast_operator
	explicit operator CSCMatrix<T>() const {
		auto n = getOrder();
		// generate pattern
		DOFsConnectivityList list(n);
		Index i, j, k;
		for (j = 0; j < n; ++j) { // for each column
			list[j].insert(j); // add diag
			for (k = _colptr[j]; k < _colptr[j + 1]; ++k) { // for each row
				i = _rowind[k];
				list[j].insert(i);
				list[i].insert(j);
			}
		}
		CSCMatrix<T> A(n, n);
		A.generatePatternFrom(list);
		// fill in values
		for (j = 0; j < n; ++j) {
			A(j, j) = _diag[j];
			for (k = _colptr[j]; k < _colptr[j + 1]; ++k) {
				i = _rowind[k];
				A(i, j) = _lval[k];
				A(j, i) = _uval[k];
			}
		}
		return A;
	}
};

// implementation

template <typename T>
CSlCMatrix<T>::CSlCMatrix(Index n, Index nnz)
	: AbstractMatrix<T>(n)
	, _colptr(n + 1)
	, _rowind(nnz)
	, _lval(nnz)
	, _uval(nnz)
	, _diag(n) {
	_colptr[n] = nnz;
}

// private methods

template <typename T>
T& CSlCMatrix<T>::_set(Index i, Index j) {
	if (i == j) return _diag[i];
	Index li = i, lj = j;
	if (i < j) std::swap(li, lj);
		for (Index k = _colptr[lj]; k < _colptr[lj + 1]; ++k)
			if (_rowind[k] == li) return i > j ? _lval[k] : _uval[k];
	throw std::invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
}

template <typename T>
T CSlCMatrix<T>::_get(Index i, Index j) const {
	if (i == j) return _diag[i];
	Index li = i, lj = j;
	if (i < j) std::swap(li, lj);
	for (Index k = _colptr[lj]; k < _colptr[lj + 1]; ++k)
		if (_rowind[k] == li) return i > j ? _lval[k] : _uval[k];
	return 0.;
}

// public methods

template <typename T>
CSlCMatrix<T>& CSlCMatrix<T>::operator=(T const & val) {
	std::fill(_lval.begin(), _lval.end(), val);
	std::fill(_uval.begin(), _uval.end(), val);
	std::fill(_diag.begin(), _diag.end(), val);
	return *this;
}

template <typename T>
void CSlCMatrix<T>::mult(T const * by, T* result) const {
	for (Index j = 0; j < _w; ++j) {
		result[j] += _diag[j] * by[j];
		for (Index i = _colptr[j]; i < _colptr[j + 1]; ++i) {
			result[_rowind[i]] += _lval[i] * by[j];
			result[j]          += _uval[i] * by[_rowind[i]];
		}
	}
}

template <typename T>
CSlCMatrix<T>& CSlCMatrix<T>::importSparse(std::istream& input) {
	// stdin structure:
	// (1) n =: _h = _w (size of square matrix)
	// (2) numb of nonzeros in lower (upper) triangular part =: _colptr[n + 1]
	// (3) n + 1    elements of _colptr,
	// (4) _colptr[n] elements of _rowind,
	// (5) _colptr[n] elements of _lval,
	// (6) _colptr[n] elements of _uval, and
	// (7) n        elements of _diag
	Index nnz;
	input >> _h >> nnz;
	_w = _h;
	_colptr.resize(_w + 1);
	_rowind.resize(nnz);
	_lval.resize(nnz);
	_uval.resize(nnz);
	_diag.resize(_w);
	input >> _colptr >> _rowind >> _lval >> _uval >> _diag;
	return *this;
}

template <typename T>
void CSlCMatrix<T>::exportSparse(std::ostream& output) const {
	output << _w << ' ' << _colptr[_w] << '\n'
	       << _colptr << '\n'
	       << _rowind << '\n'
	       << _lval << '\n'
	       << _uval << '\n'
	       << _diag;
}

template <typename T>
CSlCMatrix<T>& CSlCMatrix<T>::generatePatternFrom(DOFsConnectivityList const & list) {
	// (1) compute column pointers
	for (Index i = 0; i < _w; ++i) _colptr[i + 1] = _colptr[i] + list[i].size();
	// (2) compute row indicies
	_rowind.reserve(_colptr[_w]);
	_lval.resize(_colptr[_w]);
	_uval.resize(_colptr[_w]);
	for (Index col = 0; col < _w; ++col)
		for (auto row : list[col])
			_rowind.emplace_back(row);
	return *this;
}

template <typename T>
CSlCMatrix<T>& CSlCMatrix<T>::enforceDirichletBCs(Index2Value<T> const & ind2val, T* rhs) {
	// zero out row and col, set diag = 1, fix rhs
	decltype(ind2val.begin()) kvpIter;
	Index i, j, k;
	for (j = 0; j < getOrder(); ++j) {
		if ((kvpIter = ind2val.find(j)) != ind2val.end()) {
			for (k = _colptr[j]; k < _colptr[j + 1]; ++k) {
				i = _rowind[k];
				rhs[i] -= _lval[k] * kvpIter->second;
				_lval[k] = 0.;
				_uval[k] = 0.;
			}
			_diag[j] = 1.;
			rhs[j] = kvpIter->second;
		}
		else 
			for (k = _colptr[j]; k < _colptr[j + 1]; ++k) {
				i = _rowind[k];
				if ((kvpIter = ind2val.find(i)) != ind2val.end()) {
					rhs[j] -= _uval[k] * kvpIter->second;
					_lval[k] = 0.;
					_uval[k] = 0.;
				}
			}
	}
	return *this;
}

// find y := [ L + wD ]^-1 . x
template <typename T>
std::vector<T> CSlCMatrix<T>::forwSubst(std::vector<T> const & x, double w = 1.) const {
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
std::vector<T> CSlCMatrix<T>::backSubst(std::vector<T> const & x, double w = 1.) const {
	std::vector<T> y { x };
	SignedIndex i;
	Index j, k;
	for (i = getOrder() - 1; i >= 0; --i) {
		for (k = _colptr[i]; k < _colptr[i + 1]; ++k) {
			j = _rowind[k];
			auto& u_ij = _uval[k]; // = l_ji
			y[i] -= u_ij * y[j];
		}
		if (w) y[i] /= (w * _diag[i]);
	}
	return y;
}

template <typename T>
std::vector<T> CSlCMatrix<T>::diagSubst(std::vector<T> const & x) const {
	std::vector<T> y(getOrder());
	for (Index i = 0; i < getOrder(); ++i)
		y[i] = x[i] / _diag[i];
	return y;
}

template <typename T>
std::vector<T> CSlCMatrix<T>::multDiag(std::vector<T> const & x) const {
	std::vector<T> y(getOrder());
	for (Index i = 0; i < getOrder(); ++i)
		y[i] = x[i] * _diag[i];
	return y;
}

// ILDU(0) Crout decomposition from Saad, p. 333
// A ~ (L + I) . D . (I + U) 
template <typename T>
CSlCMatrix<T>& CSlCMatrix<T>::decompose() {
	Index n = getOrder();
	auto colptr_first = _colptr;
	std::vector<std::list<Index>> nnz_rows_of_col(n);
	std::vector<T> kth_col(n), kth_row(n);
	for (Index k = 0; k < n; ++k) {
		// copy column and row to be updated on kth step
		for (Index m = _colptr[k]; m < _colptr[k + 1]; ++m) {
			Index i = _rowind[m];
			kth_col[i] = _lval[m];
			kth_row[i] = _uval[m];
		}
		// loop over only those cols which get multiplied by a nnz row
		for (Index j : nnz_rows_of_col[k]) {
			auto u_jk = _get(j, k), l_kj = _get(k, j);
			_diag[k] -= _diag[j] * l_kj * u_jk;
			while (_rowind[colptr_first[j]] <= k) ++colptr_first[j];
			for (Index m = colptr_first[j]; m < _colptr[j + 1]; ++m) {
				Index i = _rowind[m]; auto l_ij = _lval[m], u_ji = _uval[m];
				kth_col[i] -= _diag[j] * u_jk * l_ij;
				kth_row[i] -= _diag[j] * l_kj * u_ji;
			}
		}
		// update nnz_rows_of_col lists, kth column, and kth row of the matrix
		for (Index m = _colptr[k]; m < _colptr[k + 1]; ++m) {
			Index i = _rowind[m];
			nnz_rows_of_col[i].push_back(k);
			_lval[m] = kth_col[i] / _diag[k];
			_uval[m] = kth_row[i] / _diag[k];
		}
	}
	//Index n = getOrder();
	//auto colptr_first = _colptr;
	//std::list<Index> nnz_rows;
	//std::vector<T> kth_col(n), kth_row(n);
	//for (Index k = 0; k < n; ++k) {
	//	// copy column and row to be updated on kth step
	//	for (Index m = colptr[k]; m < _colptr[k + 1]; ++m) {
	//		Index i = _rowind[m];
	//		kth_col[i] = _lval[m];
	//		kth_row[i] = _uval[m];
	//	}
	//	// loop over only those cols which get multiplied by a nnz row
	//	for (Index j : nnz_rows) { 
	//		auto u_jk = _get(j, k), l_kj = _get(k, j);
	//		_diag[k] -= _diag[j] * l_kj * u_jk;
	//		for (m = colptr_first[j]; m < _colptr[j + 1]; ++m) {
	//			Index i = _rowind[m];
	//			auto& l_ij = _lval[m];
	//			auto& u_ji = _uval[m];
	//			kth_col[i] -= _diag[j] * u_jk * l_ij;
	//			kth_row[i] -= _diag[j] * l_kj * u_ji;
	//		}
	//	}
	//	// update column and row of the matrix
	//	for (Index m = colptr[k]; m < _colptr[k + 1]; ++m) {
	//		Index i = _rowind[m];
	//		_lval[m] = kth_col[i] / _diag[k];
	//		_uval[m] = kth_row[i] / _diag[k];
	//	}
	//	// update nnz_rows and colptr_first for iteration k + 1
	//	nnz_rows.clear();
	//	for (Index j = 0; j <= k; ++j) 
	//		if (_rowind[colptr_first[j]] == k + 1) {
	//			++colptr_first[j];
	//			nnz_rows.insert(j);
	//		}
	//}
	return *this;
}