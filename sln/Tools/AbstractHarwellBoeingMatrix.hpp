#pragma once
#include "AbstractSparseMatrix.hpp"
#include "HarwellBoeingHeader.hpp"
#include "Parameters.hpp"

/*
	Alexander Žilyakov, Sep 2016
*/

template <typename T>
class AbstractHarwellBoeingMatrix
	: virtual public AbstractSparseMatrix<T> {
public:
	// construct matrix from HB–file
	virtual HarwellBoeingHeader importHarwellBoeing(std::string const &) = 0;
	// save matrix to HB–file
	virtual void exportHarwellBoeing(std::string const &, Parameters const & params = {}) const = 0; 
};		   