#pragma once
#include <algorithm> // for_each
#include <type_traits> // is_same 
#include "AbstractFEMatrix.hpp"
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
	: public AbstractFEMatrix<T>
	, public AbstractHarwellBoeingMatrix<T>
	, public AbstractTransposeMultipliableMatrix<T> {
	// …it’s similar to CSR (check out CSRMatrix.hpp for details) except 
	// _values vector contains nozero matrix values col by col, not row by row
	// so we omit details here
	std::vector<Index> _colptr, 
	                   _rowind; 
	std::vector<T>     _values;
	// virtual methods to be implemented
	T& _set(Index, Index) final;
	T  _get(Index, Index) const final;
public:
	CSCMatrix(std::vector<Index> const & colptr, std::vector<Index> const & rowind, std::vector<T> const values, Index h)
		: AbstractMatrix<T>(h, colptr.size() - 1)
		, _colptr(colptr)
		, _rowind(rowind)
		, _values(values)
	{
	}
	explicit CSCMatrix(Index h = 1, Index w = 1, Index nnz = 0);
	// virtual methods to be implemented
	Index nnz() const final { return _colptr[_w]; }
	CSCMatrix& operator=(T const & val) final;
	void mult(T const * by, T* result) const final;
	void multByTranspose(T const * by, T* result) const final;
	// i/o
	CSCMatrix& importSparse(std::istream& from = cin) final;
	void       exportSparse(std::ostream& to   = cout) const final;
	// almost same methods from base class (in order to work w/ strings instead of streams)
	using AbstractSparseMatrix::importSparse;
	using AbstractSparseMatrix::exportSparse;
	HarwellBoeingHeader importHarwellBoeing(std::string const &) final; // load from Harwell–Boeing file
	void                exportHarwellBoeing(std::string const &, Parameters const & params = {}) const final;
	CSCMatrix& generatePatternFrom(DOFsConnectivityList const &) final;
	CSCMatrix& enforceDirichletBCs(Index2Value<T> const &, T*) final;
	CSCMatrix& modifyColumn(Index j, std::function<void(T&)> const & f) {
		for (Index i = _colptr[j]; i < _colptr[j + 1]; ++i) f(_values[i]);
	}
};

// implementation

template <typename T>
CSCMatrix<T>::CSCMatrix(Index h = 1, Index w = 1, Index nnz = 0)
	: AbstractMatrix<T>(h, w)
	, _colptr(w + 1)
	, _rowind(nnz)
	, _values(nnz) {
	_colptr[w] = nnz;
}

// private methods

template <typename T>
T& CSCMatrix<T>::_set(Index i, Index j) {
	for (Index k = _colptr[j]; k < _colptr[j + 1]; ++k)
		if (_rowind[k] == i) return _values[k];
	throw std::invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
}

template <typename T>
T CSCMatrix<T>::_get(Index i, Index j) const {
	for (Index k = _colptr[j]; k < _colptr[j + 1]; ++k)
		if (_rowind[k] == i) return _values[k];
	return 0.;
}

// public methods

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::operator=(T const & val) {
	fill(_values.begin(), _values.end(), val);
	return *this;
}

template <typename T>
void CSCMatrix<T>::mult(T const * by, T* result) const {
	for (Index j = 0; j < _w; ++j)
		for (Index i = _colptr[j]; i < _colptr[j + 1]; ++i)
			result[_rowind[i]] += _values[i] * by[j];
}

template <typename T>
void CSCMatrix<T>::multByTranspose(T const * by, T* result) const {
	for (Index j = 0; j < _w; ++j) 
		for (Index i = _colptr[j]; i < _colptr[j + 1]; ++i)
			result[j] += _values[i] * by[_rowind[i]];
}

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::importSparse(std::istream& input) {
	// stdin structure:
	// (1) _h
	// (2) _w
	// (3) numb of nonzeros =: _colptr[_w]
	// (4) _w + 1    elements of _colptr,
	// (5) _colptr[_w] elements of _rowind, and
	// (6) _colptr[_w] elements of _values
	Index nnz;
	input >> _h >> _w >> nnz;
	_colptr.resize(_w + 1);
	_rowind.resize(nnz);
	_values.resize(nnz);
	input >> _colptr >> _rowind >> _values;
	return *this;
}

template <typename T>
void CSCMatrix<T>::exportSparse(std::ostream& output) const {
	output << _h << ' ' << _w << ' ' << _colptr[_w] << '\n'
	       << _colptr << '\n'
	       << _rowind << '\n'
	       << _values;
}

template <typename T>
HarwellBoeingHeader CSCMatrix<T>::importHarwellBoeing(std::string const & fileName) {
	// load header
	HarwellBoeingHeader header;
	loadHarwellBoeingHeader_f90(fileName.c_str(), &header);
	// check matrix type
	if (header.mxtype[1] != 'U' && header.mxtype[1] != 'R') throw std::invalid_argument("CSC type is useful for rect / unsymmetric pattern (.*r*, .*u*) Harwell-Boeing matrices; try SymmetricCSlC for (.*s*) type");
	if (header.mxtype[2] != 'A') throw std::invalid_argument("only assembled (.**a) Harwell-Boeing matrices are supported");
	if (std::is_same<T, double>::value) { // real matrix 
		if (header.mxtype[0] == 'C') throw std::invalid_argument("instantiate matrix w/ complex<double> in order to work w/ complex (.c**) Harwell-Boeing matrices");
	}
	else // complex matrix
		if (header.mxtype[0] == 'R') throw std::invalid_argument("instantiate matrix w/ double in order to work w/ real (.r**) Harwell-Boeing matrices");
	// fix internal structure
	_h = header.nrow;
	_w = header.ncol;
	_colptr.resize(_w + 1);
	_rowind.resize(header.nnzero);
	_values.resize(header.nnzero);
	// load structure and values
	loadHarwellBoeingStruct_f90(fileName.c_str(), &header, _colptr.data(), _rowind.data(), 
	                            reinterpret_cast<double*>(_values.data())); // http://en.cppreference.com/w/cpp/numeric/complex > non-static data members
	// TODO: FORTRAN should care for this
	for_each(_colptr.begin(), _colptr.end(), [](Index& i) { --i; });
	for_each(_rowind.begin(), _rowind.end(), [](Index& i) { --i; });

	return header;
}

template <typename T>
void CSCMatrix<T>::exportHarwellBoeing(std::string const & fileName, Parameters const & params) const {
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
		"CRA", _h, _w, _colptr[_w], 0 /* always zero for assembled matrices */,
		// line 4
		"", "", "(3E26.16)", ""
		// w/o line 5 — no RHS
	};
	// fix line 1
	if (params.find("title") != params.end()) strcpy_s(header.title, params.at("title").c_str());
	if (params.find("key")   != params.end()) strcpy_s(header.key,   params.at("key").c_str());
	// fix lines 2—4
	if (_w == _h) header.mxtype[1] = 'U'; // unsymmetric matrix
	auto numbOfDigits = [](Index n) { // helper
		return floor(log10(n ? n : 1) + 1);
	};
	// fix format of column pointers 
	Index numbOfRecs = _w + 1,	// numb of records for col pointers (row indicies, values)
	      recLength = numbOfDigits(_colptr[_w] + 1 /* = max value of _colptr is its last element */),	// length of single record
	      recsPerLine = 80 /* = numb of cols in HB file */ / recLength; // numb of records per one line
	header.ptrcrd = ceil( (float) numbOfRecs / recsPerLine); // numb of lines for column pointers
	auto fmt = '(' + std::to_string(recsPerLine) + "I" + std::to_string(recLength) + ')'; // FORTRAN format
	strcpy_s(header.ptrfmt, fmt.c_str());
	// fix format of row indicies
	numbOfRecs = _colptr[_w];
	recLength = numbOfDigits(_h /* = max value of _rowind is numb of rows */);
	recsPerLine = 80 / recLength;
	header.indcrd = ceil( (float) numbOfRecs / recsPerLine ); // numb of lines for row indicies
	fmt = '(' + std::to_string(recsPerLine) + "I" + std::to_string(recLength) + ')';
	strcpy_s(header.indfmt, fmt.c_str());
	if (params.find("pattern") != params.end()) { // pattern only
		header.mxtype[0] = 'P';
		header.valfmt[0] = '\0';
	}
	else if (std::is_same<T, double>::value) { // real matrix
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
	saveHarwellBoeingStruct_f90(fileName.c_str(), &header, _colptr.data(), _rowind.data(), reinterpret_cast<double const *>(_values.data()));
}

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::generatePatternFrom(DOFsConnectivityList const & list) {
	// (1) compute column pointers
	for (Index i = 0; i < _w; ++i) _colptr[i + 1] = _colptr[i] + list[i].size();
	// (2) compute row indicies
	_rowind.reserve(_colptr[_w]);
	_values.resize(_colptr[_w]);
	for (Index col = 0; col < _w; ++col)
		for (auto row : list[col])
			_rowind.emplace_back(row);
	return *this;
}

template <typename T>
CSCMatrix<T>& CSCMatrix<T>::enforceDirichletBCs(Index2Value<T> const & ind2val, T* rhs) {
	// zero out column and fix rhs
	// wrn: works for non-square matrices, we do not set diag values to unities
	// if you want so, use CSlC / symmetric CSlC matrix instead
	decltype(ind2val.begin()) kvpIter;
	Index i, j, k;
	for (j = 0; j < _w; ++j) 
		if ((kvpIter = ind2val.find(j)) != ind2val.end())
			for (k = _colptr[j]; k < _colptr[j + 1]; ++k) {
				i = _rowind[k];
				rhs[i] -= _values[k] * kvpIter->second;
				_values[k] = 0.;
			}
	return *this;
}