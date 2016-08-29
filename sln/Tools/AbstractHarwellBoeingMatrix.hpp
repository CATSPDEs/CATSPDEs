#pragma once
#include "AbstractSparseMatrix.hpp"
#include "HarwellBoeingHeader.hpp"

template <typename T>
class AbstractHarwellBoeingMatrix
	: public AbstractSparseMatrix<T> {
protected:
	HarwellBoeingHeader* _header;
public:
	AbstractHarwellBoeingMatrix() : _header(nullptr) {}
	AbstractHarwellBoeingMatrix(HarwellBoeingHeader& header) : _header(&header) {}
	virtual void loadHarwellBoeing(string const &) = 0; // construct matrix from HB–file
	virtual void saveHarwellBoeing() const = 0; // save matrix to HB–file
};		   

// these routines will be linked from
// Harwell–Boeing IO (FORTRAN Static Library) project
// so you have to include .lib file from the above project
// in order to use CSC–matrix in your code
//extern "C" {
//	// read header
//	void readHarwellBoeingHeader(char const *, HarwellBoeingHeader*);
//	// read struct and values
//	void readHarwellBoeingStruct(char const *, HarwellBoeingHeader const *, size_t*, size_t*, double*); 
//	// ″                     , but for symmetric matrices
//	// HB stores diagonal in lower triangle, but we in CATSPDEs prefer store it
//	// separately since in FEM / FVM diags are nonzero
//	// so we need slightly modified version of the above func for symmetric matrices
//	void readHarwellBoeingStructSym(char const *, HarwellBoeingHeader const *, size_t*, size_t*, double*, double*);
//	void writeHarwellBoeing     (char const *, HarwellBoeingHeader const *, size_t const *, size_t const *, double const *);
//}