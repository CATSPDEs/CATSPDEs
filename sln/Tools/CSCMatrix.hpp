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
	CSCMatrix(HarwellBoeingHeader&);
	// virtual methods to be implemented
	CSCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	istream& loadSparse(istream& from = cin) final;
	ostream& saveSparse(ostream& to   = cout) const final;
	void loadHarwellBoeing(string const &) final; // load from Harwell–Boeing file
	void saveHarwellBoeing() const final;
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
CSCMatrix<T>::CSCMatrix(HarwellBoeingHeader& header)
	: AbstractMatrix<T>(header.nrow, header.ncol)
	, AbstractHarwellBoeingMatrix<T>(header)
	, _jptr(header.ncol + 1)
	, _iptr(header.nnzero)
	, _mval(header.nnzero) {
	_jptr[header.ncol] = header.nnzero;
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

//template <typename T>
//void CSCMatrix<T>::logStructure(SingletonLogger& logger, unsigned short min2print) const {
//	//// _jptr
//	//size_t n = _w + 1; // actual numb of elements
//	//unsigned short m = min(min2print, n); // numb of elements to print
//	//logger.log("_jptr, vector of column pointers");
//	//for (i = 0; i < m; ++i) logger.mes("[" + to_string(i) + "] -> " + to_string(_jptr[i]));
//	//if (m < n) {
//	//	if (m < n - 1) logger.mes(". . .");
//	//	logger.mes("[" + to_string(n - 1) + "] -> " + to_string(_jptr[n - 1]));
//	//}
//	//// _iptr
//	//n = _jptr[_w];
//	//m = min(min2print, n);
//	//logger.log("_iptr, vector of row indicies");
//	//for (i = 0; i < m; ++i) logger.mes("[" + to_string(i) + "] -> " + to_string(_iptr[i]));
//	//if (m < n) {
//	//	if (m < n - 1) logger.mes(". . .");
//	//	logger.mes("[" + (n - 1) + "] -> " + to_string(_iptr[n - 1]));
//	//}
//	//// _mval
//	//auto to_string = [](complex<double> const & c) { // helper for complex T
//	//	return std::to_string(c.real()) + " + i * " + std::to_string(c.imag());
//	//};
//	//logger.log("_mval, vector matrix values");
//	//for (i = 0; i < m; ++i) logger.mes("[" + to_string(i) + "] -> " + to_string(_mval[i]));
//	//if (m < n) {
//	//	if (m < n - 1) logger.mes(". . .");
//	//	logger.mes("[" + (n - 1) + "] -> " + to_string(_mval[n - 1]));
//	//}
//	//if (is_same<T, double>::value) { // real
//	//	for (i = 0; i < n; ++i) logger.mes("mval[" + to_string(i) + "] = " + to_string(mval[i]));
//	//	logger.mes("mval[" + to_string(header.nnzero - 1) + "] = " + to_string(mval[header.nnzero - 1]));
//	//}
//	//else { // complex
//	//	for (i = 0; i < n; ++i) logger.mes("mval[" + to_string(i) + "] = " + to_string(mval[2 * i]) + " + i * " + to_string(mval[2 * i + 1]));
//	//	logger.mes("mval[" + to_string(header.nnzero - 1) + "] = " + to_string(mval[2 * header.nnzero - 2]) + " + i * " + to_string(mval[2 * header.nnzero - 1]));
//	//}
//}

extern "C" void loadHarwellBoeingStruct_f90(char const *, HarwellBoeingHeader const *, size_t*, size_t*, double*);

template <typename T>
void CSCMatrix<T>::loadHarwellBoeing(string const & fileName) {
	// check matrix type
	if (_header->mxtype[1] != 'U' && _header->mxtype[1] != 'R') throw invalid_argument("CSC type is useful for rect / unsymmetric pattern (.*r*, .*u*) Harwell-Boeing matrices; try SymmetricCSlR");
	if (_header->mxtype[2] != 'A') throw invalid_argument("only assembled (.**a) Harwell-Boeing matrices are supported");
	if (is_same<T, double>::value) { // real matrix 
		if (_header->mxtype[0] == 'C') throw invalid_argument("instantiate matrix w/ complex<double> in order to work w/ complex (.c**) Harwell-Boeing matrices");
	}
	else // complex matrix
		if (_header->mxtype[0] == 'R') throw invalid_argument("instantiate matrix w/ double in order to work w/ real (.r**) Harwell-Boeing matrices");
	// read structure and values
	double* mval = reinterpret_cast<double*>(_mval.data()); // http://en.cppreference.com/w/cpp/numeric/complex > non-static data members
	loadHarwellBoeingStruct_f90(fileName.c_str(), _header, _jptr.data(), _iptr.data(), mval);
	for_each(_jptr.begin(), _jptr.end(), [](size_t& i) { --i; });
	for_each(_iptr.begin(), _iptr.end(), [](size_t& i) { --i; });
}

template <typename T>
void CSCMatrix<T>::saveHarwellBoeing(/*string const & fileName, string const & title, string const & key*/) const {
	// …
}