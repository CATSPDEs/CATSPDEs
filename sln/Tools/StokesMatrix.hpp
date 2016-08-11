#pragma once
#include "SymmetricCSlRMatrix.hpp"
#include "CSRMatrix.hpp"

class StokesMatrix 
	: public AbstractSparseMatrix<double>
	, public IMultipliable<double> {
	SymmetricCSlRMatrix<double> _A11, _A22;
	CSRMatrix<double> _B1, _B2;

	vector<double> _mult(vector<double> const &) override;

public:
};