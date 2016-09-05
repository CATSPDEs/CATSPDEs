/*
	Alexander Žilyakov, Sep 2016
*/
#pragma once
#include "AbstractSparseMatrix.hpp"
#include "HarwellBoeingHeader.hpp"
#include "Parameters.hpp"

template <typename T>
class AbstractHarwellBoeingMatrix
	: public AbstractSparseMatrix<T> {
public:
	virtual AbstractHarwellBoeingMatrix& loadHarwellBoeing(HarwellBoeingHeader const &, string const &) = 0; // construct matrix from HB–file
	virtual void                         saveHarwellBoeing(string const &, Parameters const & params = {}) const = 0; // save matrix to HB–file
};		   