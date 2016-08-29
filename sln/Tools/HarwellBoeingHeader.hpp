#pragma once
#include <iostream>

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

std::ostream& operator<<(std::ostream& output, HarwellBoeingHeader const & header) {
	output << "line 1\n"
		<< "\ttitle:  " << header.title << '\n'
		<< "\tkey:    " << header.key << '\n'
		<< "line 2\n"
		<< "\ttotcrd: " << header.totcrd << '\n'
		<< "\tptrcrd: " << header.ptrcrd << '\n'
		<< "\tindcrd: " << header.indcrd << '\n'
		<< "\tvalcrd: " << header.valcrd << '\n'
		<< "\trhscrd: " << header.rhscrd << '\n'
		<< "line 3\n"
		<< "\tmxtype: " << header.mxtype << '\n'
		<< "\tnrow:   " << header.nrow << '\n'
		<< "\tncol:   " << header.ncol << '\n'
		<< "\tnnzero: " << header.nnzero << '\n'
		<< "\tneltvl: " << header.neltvl << '\n'
		<< "line 4\n"
		<< "\tptrfmt: " << header.ptrfmt << '\n'
		<< "\tindfmt: " << header.indfmt << '\n'
		<< "\tvalfmt: " << header.valfmt << '\n'
		<< "\trhsfmt: " << header.rhsfmt << '\n';
	if (header.rhscrd > 0)
		output << "line 5\n"
		<< "\trhstyp: " << header.rhstyp << '\n'
		<< "\tnrhs: " << header.nrhs << '\n'
		<< "\tnrhsix: " << header.nrhsix << '\n';
	return output;
}

// these routines will be linked from
// Harwell–Boeing IO (FORTRAN Static Library) project
// so you have to include .lib file from the above project
// in order to use CSC–like matrices in your code
extern "C" {
	void loadHarwellBoeingHeader_f90(char const *, HarwellBoeingHeader*);
	void saveHarwellBoeingHeader_f90();
}