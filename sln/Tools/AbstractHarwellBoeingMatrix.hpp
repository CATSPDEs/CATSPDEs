#pragma once
#include "AbstractSparseMatrix.hpp"
#include "SingletonLogger.hpp"

struct HarwellBoeingHeader;

template <typename T>
class AbstractHarwellBoeingMatrix
	: public AbstractSparseMatrix<T> {
public:
	// we should also provide
	// AbstractHarwellBoeingMatrixDERIVED(string const &, SingletonLogger&);
	// in derived classe
	virtual void saveHarwellBoeing(string const &, string const &, string const &) = 0;
};

// interoperable C —> FORTRAN structure 
struct HarwellBoeingHeader {
	// line 1
	char title[73], // +1 element for null–character of C–strings (i.e. 73 instead of 72 needed,
	     key[9];    // 9 instead of 8 etc.) 
	// line 2
	size_t totcrd, // total numb of lines excluding header
	       ptrcrd, //       ″             for col pointers
	       indcrd, //       ″                 row indicies
	       valcrd, //       ″                 numerical values of matrix
	       rhscrd; //       ″                 right–hand sides
	// line 3
	char mxtype[4]; // matrix type
	size_t nrow,    // numb of rows
           ncol,    // ″       cols
           nnzero,  // ″       row indicies / entries in assembled matrix
           neltvl;  // ″       elemental matrix entries
	// line 4
	char ptrfmt[17], // FORTRAN format for col pointers
	     indfmt[17], // ″                  row indicies
	     valfmt[21], // ″                  numerical values of matrix
	     rhsfmt[21]; // ″                  values right–hand sides
	// line 5 (optional)
	char rhstyp[4]; // right–hand side type
	size_t nrhs,    // numb of right–hand sides
	       nrhsix;  // ″       row indicies
};

void logHarwellBoeingHeader(HarwellBoeingHeader const & header, SingletonLogger& logger) {
	logger.log("line 1");
	logger.mes("title:  " + string(header.title));
	logger.mes("key:    " + string(header.key));
	logger.log("line 2");
	logger.mes("totcrd: " + to_string(header.totcrd));
	logger.mes("ptrcrd: " + to_string(header.ptrcrd));
	logger.mes("indcrd: " + to_string(header.indcrd));
	logger.mes("valcrd: " + to_string(header.valcrd));
	logger.mes("rhscrd: " + to_string(header.rhscrd));
	logger.log("line 3");
	logger.mes("mxtype: " + string(header.mxtype));
	logger.mes("nrow:   " + to_string(header.nrow));
	logger.mes("ncol:   " + to_string(header.ncol));
	logger.mes("nnzero: " + to_string(header.nnzero));
	logger.mes("neltvl: " + to_string(header.neltvl));
	logger.log("line 4");
	logger.mes("ptrfmt: " + string(header.ptrfmt));
	logger.mes("indfmt: " + string(header.indfmt));
	logger.mes("valfmt: " + string(header.valfmt));
	logger.mes("rhsfmt: " + string(header.rhsfmt));
	if (header.rhscrd > 0) {
		logger.log("line 5");
		logger.mes("rhstyp: " + string(header.rhstyp));
		logger.mes("nrhs: " + to_string(header.nrhs));
		logger.mes("nrhsix: " + to_string(header.nrhsix));
	}
}

// these routines will be linked from
// Harwell–Boeing IO (FORTRAN Static Library) project
// so you have to include .lib file from the above project
// in order to use CSC–matrix in your code
extern "C" {
	void readHarwellBoeingHeader(char const *, HarwellBoeingHeader*);
	void readHarwellBoeingStruct(char const *, HarwellBoeingHeader const *, size_t*, size_t*, double*);
	void writeHarwellBoeing     (char const *, HarwellBoeingHeader const *, size_t const *, size_t const *, double const *);
}