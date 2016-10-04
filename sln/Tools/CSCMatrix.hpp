#pragma once
#include <algorithm> // for_each
#include <type_traits> // is_same 
#include "AbstractHarwellBoeingMatrix.hpp"
#include "AbstractTransposeMultipliableMatrix.hpp"

/*
	Alexander Žilyakov, Sep 2016
*/

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
	vector<size_t> _jptr, 
	               _iptr; 
	vector<T>      _mval;
	// virtual methods to be implemented
	size_t _nnz() const final { return _jptr[_w]; }
	T& _set(size_t, size_t) final;
	T  _get(size_t, size_t) const final;
public:
	explicit CSCMatrix(size_t h = 1, size_t w = 1, size_t nnz = 0);
	~CSCMatrix() {}
	// virtual methods to be implemented
	CSCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	CSCMatrix& loadSparse(istream& from = cin) final;
	void       saveSparse(ostream& to   = cout) const final;
	CSCMatrix& loadHarwellBoeing(string const &, HarwellBoeingHeader* headerPtr = nullptr) final; // load from Harwell–Boeing file
	void       saveHarwellBoeing(string const &, Parameters const & params = {}) const final;
};

// implementation

template <typename T>
CSCMatrix<T>::CSCMatrix(size_t h = 1, size_t w = 1, size_t nnz = 0)
	: AbstractMatrix<T>(h, w)
	, _jptr(w + 1)
	, _iptr(nnz)
	, _mval(nnz) {
	_jptr[w] = nnz;
}

// private methods

template <typename T>
T& CSCMatrix<T>::_set(size_t i, size_t j) {
	for (size_t k = _jptr[j]; k < _jptr[j + 1]; ++k)
		if (_iptr[k] == i) return _mval[k];
		else if (_iptr[k] > i) throw invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
	throw invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
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
CSCMatrix<T>& CSCMatrix<T>::loadSparse(istream& input) {
	// stdin structure:
	// (1) _h
	// (2) _w
	// (3) numb of nonzeros =: _jptr[_w]
	// (4) _w + 1    elements of _jptr,
	// (5) _jptr[_w] elements of _iptr, and
	// (6) _jptr[_w] elements of _mval
	size_t nnz;
	input >> _h >> _w >> nnz;
	_jptr.resize(_w + 1);
	_iptr.resize(nnz);
	_mval.resize(nnz);
	input >> _jptr >> _iptr >> _mval;
	return *this;
}

template <typename T>
void CSCMatrix<T>::saveSparse(ostream& output) const {
	output << _h << ' ' << _w << ' ' << _jptr[_w] << '\n'
	       << _jptr << '\n'
	       << _iptr << '\n'
	       << _mval;
}

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::loadHarwellBoeing(string const & fileName, HarwellBoeingHeader* headerPtr) {
	// load header
	HarwellBoeingHeader header;
	loadHarwellBoeingHeader_f90(fileName.c_str(), &header);
	// check matrix type
	if (header.mxtype[1] != 'U' && header.mxtype[1] != 'R') throw invalid_argument("CSC type is useful for rect / unsymmetric pattern (.*r*, .*u*) Harwell-Boeing matrices; try SymmetricCSlC for (.*s*) type");
	if (header.mxtype[2] != 'A') throw invalid_argument("only assembled (.**a) Harwell-Boeing matrices are supported");
	if (is_same<T, double>::value) { // real matrix 
		if (header.mxtype[0] == 'C') throw invalid_argument("instantiate matrix w/ complex<double> in order to work w/ complex (.c**) Harwell-Boeing matrices");
	}
	else // complex matrix
		if (header.mxtype[0] == 'R') throw invalid_argument("instantiate matrix w/ double in order to work w/ real (.r**) Harwell-Boeing matrices");
	// fix internal structure
	_h = header.nrow;
	_w = header.ncol;
	_jptr.resize(_w + 1);
	_iptr.resize(header.nnzero);
	_mval.resize(header.nnzero);
	// load structure and values
	loadHarwellBoeingStruct_f90(fileName.c_str(), &header, _jptr.data(), _iptr.data(), 
	                            reinterpret_cast<double*>(_mval.data())); // http://en.cppreference.com/w/cpp/numeric/complex > non-static data members
	// TODO: FORTRAN should care for this
	for_each(_jptr.begin(), _jptr.end(), [](size_t& i) { --i; });
	for_each(_iptr.begin(), _iptr.end(), [](size_t& i) { --i; });
	// optionally return header
	if (headerPtr) *headerPtr = header;
	return *this;
}

template <typename T>
void CSCMatrix<T>::saveHarwellBoeing(string const & fileName, Parameters const & params) const {
	/* 
	available params:
		"title"   -> "Created w/ CATSPDEs: https://github.com/CATSPDEs"  matrix title (see HB–format docs), 72 chars max
		"key"     -> "CATSPDEs"                                          matrix key, 8 chars max
		"pattern"                                                        save only pattern (.**p)
	*/
	// (1) generate header
	HarwellBoeingHeader header = {
		// line 1
		"Created w/ CATSPDEs: https://github.com/CATSPDEs", "CATSPDEs",
		// line 2
		0, 0, 0, 0, 0,
		// line 3
		"CRA", _h, _w, _jptr[_w], 0 /* always zero for assembled matrices */,
		// line 4
		"", "", "(3E26.16)", ""
		// w/o line 5 — no RHS
	};
	// fix line 1
	if (params.find("title") != params.end()) strcpy_s(header.title, params.at("title").c_str());
	if (params.find("key")   != params.end()) strcpy_s(header.key,   params.at("key").c_str());
	// fix lines 2—4
	if (_w == _h) header.mxtype[1] = 'U'; // unsymmetric matrix
	auto numbOfDigits = [](size_t n) { // helper
		return floor(log10(n ? n : 1) + 1);
	};
	// fix format of column pointers 
	size_t numbOfRecs = _w + 1,	// numb of records for col pointers (row indicies, values)
	       recLength = numbOfDigits(_jptr[_w] + 1 /* = max value of _jptr is its last element */),	// length of single record
	       recsPerLine = 80 /* = numb of cols in HB file */ / recLength; // numb of records per one line
	header.ptrcrd = ceil( (float) numbOfRecs / recsPerLine); // numb of lines for column pointers
	string fmt = '(' + to_string(recsPerLine) + "I" + to_string(recLength) + ')'; // FORTRAN format
	strcpy_s(header.ptrfmt, fmt.c_str());
	// fix format of row indicies
	numbOfRecs = _jptr[_w];
	recLength = numbOfDigits(_h /* = max value of _iptr is numb of rows */);
	recsPerLine = 80 / recLength;
	header.indcrd = ceil( (float) numbOfRecs / recsPerLine ); // numb of lines for row indicies
	fmt = '(' + to_string(recsPerLine) + "I" + to_string(recLength) + ')';
	strcpy_s(header.indfmt, fmt.c_str());
	if (params.find("pattern") != params.end()) { // pattern only
		header.mxtype[0] = 'P';
		header.valfmt[0] = '\0';
	}
	else if (is_same<T, double>::value) { // real matrix
		header.mxtype[0] = 'R'; 
		// recLength   = 26 symbols — double precision: {+|-}d.dddddddddddddddde{+|-}ddddd
		// recsPerLine = 80 / 26 = 3
		header.valcrd = ceil(numbOfRecs / 3.); // numb of lines for values
	}
	else header.valcrd = ceil(2. * numbOfRecs / 3.); // for complex values
	header.totcrd = header.ptrcrd + header.indcrd + header.valcrd;
	// (2) save header
	saveHarwellBoeingHeader_f90(fileName.c_str(), &header);
	// (3) save struct
	saveHarwellBoeingStruct_f90(fileName.c_str(), &header, _jptr.data(), _iptr.data(), reinterpret_cast<double const *>(_mval.data()));
}