#pragma once
#include <algorithm> // for_each
#include <type_traits> // is_same 
#include "AbstractHarwellBoeingMatrix.hpp"
#include "AbstractTransposeMultipliableMatrix.hpp"


template <typename T> class CSCMatrix;
// matrix is stored in CSC (compressed sparse column) format…

template <typename T> using HBMatrix = CSCMatrix<T>;
// …aka Harwell–Boeing format…

template <typename T>
class CSCMatrix
	: public AbstractHarwellBoeingMatrix<T>
	, public AbstractTransposeMultipliableMatrix<T> {
	// …it’s similar to CSR (check out CSRMatrix.hpp for details) except 
	// _mval vector contains nozero matrix values col by col, not row by row
	// so we omit details here
	vector<T>      _mval; 
	vector<size_t> _jptr, 
	               _iptr; 
	// virtual methods to be implemented
	size_t _nnz() const final { return _jptr[_w]; }
	T& _set(size_t, size_t) final;
	T  _get(size_t, size_t) const final;
public:
	CSCMatrix(size_t h, size_t w, size_t nnz);
	CSCMatrix(string const &, SingletonLogger&); // load from Harwell–Boeing file
	// virtual methods to be implemented
	CSCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	istream& loadSparse(istream& from = cin) final;
	ostream& saveSparse(ostream& to   = cout) const final;
	void saveHarwellBoeing(string const &, string const &, string const &) final;
};

// implementation

template <typename T>
CSCMatrix<T>::CSCMatrix(size_t h, size_t w, size_t nnz)
	: AbstractMatrix<T>(h, w)
	, _jptr(w + 1)
	, _iptr(nnz)
	, _mval(nnz) {
	_jptr[w] = nnz;
}

template <typename T>
CSCMatrix<T>::CSCMatrix(string const & fileName, SingletonLogger& logger)
	: AbstractMatrix(1, 1) {
	HarwellBoeingHeader header;
	logger.beg("read header info from " + fileName);
		readHarwellBoeingHeader(fileName.c_str(), &header);
		logHarwellBoeingHeader(header, logger);
	logger.end();
	// check matrix type
	if (header.mxtype[1] != 'U' && header.mxtype[1] != 'R') throw invalid_argument("CSC type is useful for rect / unsymmetric pattern (.*r*, .*u*) Harwell-Boeing matrices; try SymmetricCSlR");
	if (header.mxtype[2] != 'A') throw invalid_argument("only assembles (.**a) Harwell-Boeing matrices are supported");
	// update structure
	_h = header.nrow;
	_w = header.ncol;
	_jptr.resize(header.ncol + 1);
	_iptr.resize(header.nnzero);
	_mval.resize(header.nnzero);
	logger.beg("read matrix structure / values from " + fileName);
		double* mval = reinterpret_cast<double*>(_mval.data());
		readHarwellBoeingStruct(fileName.c_str(), &header, _jptr.data(), _iptr.data(), mval);
		logger.log("decrement iptr and jptr values");
		for_each(_jptr.begin(), _jptr.end(), [](size_t& i) { --i; });
		for_each(_iptr.begin(), _iptr.end(), [](size_t& i) { --i; });
		logger.log("print internal data");
		// jptr
		short n = min(4, header.ncol + 1), i;
		for (i = 0; i < n; ++i) logger.mes("jptr[" + to_string(i) + "] = " + to_string(_jptr[i]));
		logger.mes("jptr[" + to_string(header.ncol) + "] = " + to_string(_jptr[header.ncol]));
		// iptr
		n = min(4, header.nnzero);
		for (i = 0; i < n; ++i) logger.mes("iptr[" + to_string(i) + "] = " + to_string(_iptr[i]));
		logger.mes("iptr[" + to_string(header.nnzero - 1) + "] = " + to_string(_iptr[header.nnzero - 1]));
		// mval
		if (is_same<T, double>::value) { // real
			for (i = 0; i < n; ++i) logger.mes("mval[" + to_string(i) + "] = " + to_string(mval[i]));
			logger.mes("mval[" + to_string(header.nnzero - 1) + "] = " + to_string(mval[header.nnzero - 1]));
		}
		else { // complex
			for (i = 0; i < n; ++i) logger.mes("mval[" + to_string(i) + "] = " + to_string(mval[2 * i]) + " + i * " + to_string(mval[2 * i + 1]));
			logger.mes("mval[" + to_string(header.nnzero - 1) + "] = " + to_string(mval[2 * header.nnzero - 2]) + " + i * " + to_string(mval[2 * header.nnzero - 1]));
		}
	logger.end();
}

// private methods

template <typename T>
T& CSCMatrix<T>::_set(size_t i, size_t j) {
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _mval[k];
		else if (_iptr[k] > i) throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	throw invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
}

template <typename T>
T CSCMatrix<T>::_get(size_t i, size_t j) const {
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _mval[k];
		else if (_iptr[k] > i) return 0.;
	return 0.;
}

// public methods

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::operator=(T const & val) {
	fill(_mval.begin(), _mval.end(), val);
	return *this;
}

template <typename T>
void CSCMatrix<T>::mult(T const * by, T* result) const {
	size_t i, j;
	fill(result, result + _h, 0.); // clear resulting vector
	for (j = 0; j < _w; ++j)
		for (i = _jptr[j]; i < _jptr[j + 1]; ++i)
			result[_iptr[i]] += _mval[i] * by[j];
}

template <typename T>
void CSCMatrix<T>::multByTranspose(T const * by, T* result) const {
	size_t i, j;
	for (j = 0; j < _w; ++j) {
		result[j] = 0.;
		for (i = _jptr[j]; i < _jptr[j + 1]; ++i)
			result[j] += _mval[i] * by[_iptr[i]];
	}
}

template <typename T>
istream& CSCMatrix<T>::loadSparse(istream& input) {
	// stdin structure:
	// (1) _w + 1    elements of _jptr,
	// (2) _jptr[_w] elements of _iptr, and
	// (3) _jptr[_w] elements of _mval
	return input >> _jptr >> _iptr >> _mval;
}

template <typename T>
ostream& CSCMatrix<T>::saveSparse(ostream& output) const {
	return output << _h << ' ' << _w << ' ' << _jptr[_w] << '\n'
	              << _jptr << '\n'
	              << _iptr << '\n'
	              << _mval;
}

template <typename T>
void CSCMatrix<T>::saveHarwellBoeing(string const & fileName, string const & title, string const & key) {
	// …
}