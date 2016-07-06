#include"CSlRMatrix.hpp"
CSlRMatrix::CSlRMatrix(size_t n, size_t nnz)
	: AbstractSparseMatrix(n)
	, _diag(n)
	, _iptr(n + 1)
	, _jptr(nnz)
	, _lval(nnz) 
	, _uval(nnz){
	_iptr[n] = nnz;
}

double& CSlRMatrix::_set(size_t i, size_t j) {
	if (i == j) return _diag[i];
	if (i < j) {
		for (size_t k = _iptr[j]; k < _iptr[j + 1]; ++k)
			if (_jptr[k] == i) return _uval[k];
			else if (_jptr[k] > i) throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
		throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	}
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _lval[k];
		else if (_jptr[k] > j) throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}
double CSlRMatrix::_get(size_t i, size_t j) const {
	if (i == j) return _diag[i];
	if (i < j) {
		for (size_t k = _iptr[j]; k < _iptr[j + 1]; ++k)
			if (_jptr[k] == i) return _uval[k];
			else if (_jptr[k] > i) return 0.;
		return 0.;
	}
	for (size_t k = _iptr[i]; k < _iptr[i + 1]; ++k)
		if (_jptr[k] == j) return _lval[k];
		else if (_jptr[k] > j) return 0.;
	return 0.;
}

std::istream & CSlRMatrix::loadSparse(std::istream& input) {
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
		input >> _uval[i];
	for (i = 0; i < max; ++i)
		input >> _lval[i];
	max = _n;
	for (i = 0; i < max; ++i)
		input >> _diag[i];
	return input;
}

std::ostream& CSlRMatrix::saveSparse(std::ostream& output) const {
	return output << _n << ' ' << _iptr[_n] << '\n'
		<< std::vector<size_t>(_iptr.begin(), _iptr.end() - 1) << '\n'
		// we do not want last element to be written—we already have it
		<< _jptr << '\n'
		<< _uval << '\n'
		<< _lval << '\n'
		<< _diag;
}

std::vector<double> CSlRMatrix::solve(std::vector<double> const &)
{
	//placeholder
	return std::vector<double>();
}

std::vector<double> CSlRMatrix::mult(std::vector<double> const &u) const { // compute v := A.u
	std::vector<double> v(_n, 0.);
	for (size_t i = 0; i < _n; ++i) {
		v[i] += u[i] * _diag[i];
		for (size_t j = _iptr[i]; j < _iptr[i + 1]; ++j) {
			v[i] += _lval[j] * u[_jptr[j]];
			v[_jptr[j]] += _uval[j] * u[i];
		}
	}
	return v;
}

std::vector<double> CSlRMatrix::multt(std::vector<double> const &u) const { // compute v := A.u
	std::vector<double> v(_n, 0.);
	for (size_t i = 0; i < _n; ++i) {
		v[i] += u[i] * _diag[i];
		for (size_t j = _iptr[i]; j < _iptr[i + 1]; ++j) {
			v[i] += _uval[j] * u[_jptr[j]];
			v[_jptr[j]] += _lval[j] * u[i];
		}
	}
	return v;
}

std::vector<double> operator& (CSlRMatrix const & A, std::vector<double> const & b) {
	return A.multt(b);
}

CSlRMatrix CSlRMatrix::ILU() {
	CSlRMatrix LU(*this);
	for (size_t k = 1;k < _n;k++) {
		for (size_t j = _iptr[k];j < _iptr[k + 1];j++) {
			for (size_t i = _iptr[k];i < j;i++)
				LU._lval[j] -= LU._lval[i] * 
				[&](size_t q, size_t c) {
				int p;
				for (p = _iptr[c]; p < _iptr[c + 1];p++)
					if (_jptr[p] == q) return LU._uval[p];
				return 0.;
			}(k, _jptr[j]);
			LU._lval[j] /= LU._diag[_jptr[j]];
		}
		for (size_t j = _iptr[k];j < _iptr[k + 1];j++)
			LU._diag[k] -= LU._uval[j] * LU._lval[j];
		for (size_t j = _iptr[k];j < _iptr[k + 1];j++)
			for (size_t i = _iptr[k];i < j;i++)
				LU._uval[j] -= LU._uval[i] * [&](size_t q, size_t c) {
				int p;
				for (p = _iptr[c]; p < _iptr[c + 1];p++)
					if (_jptr[p] == q) return LU._lval[p];
				return 0.;
			}(k, _jptr[j]);

	}
	return LU;
}

std::vector<double> CSlRMatrix::forwardSubstitution(std::vector<double> const &f) {
	std::vector<double> temp(f);
	std::vector<double> res(_n);
	for (int i = 0;i < _n;i++) {
		for (int j = _iptr[i];j < _iptr[i + 1];j++)
			temp[i] -= res[_jptr[j]] * _lval[j];
		res[i] = temp[i];
	}
	return res;
}

std::vector<double> CSlRMatrix::backwardSubstitution(std::vector<double> const &f) {
	std::vector<double> temp(f);
	std::vector<double> res(_n);
	for (int i = _n-1;i >= 0;i--) {
		res[i] = temp[i] / _diag[i];
		for (int j = _iptr[i];j < _iptr[i + 1];j++)
			temp[_jptr[j]] -= res[i] * _uval[j];

	}
	return res;
}