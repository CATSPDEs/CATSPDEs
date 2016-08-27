#pragma once

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

// these routines will be linked from
// Harwell–Boeing IO (FORTRAN Static Library) project
// so you have to include .lib file from the above project
// in order to use CSC–matrix in your code
extern "C" {
	void readHarwellBoeingHeader(char const *, HarwellBoeingHeader*);
	void readHarwellBoeingStruct(char const *, HarwellBoeingHeader const *, size_t*, size_t*, double*);
}