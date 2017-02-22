#pragma once
#include "AbstractMultipliableMatrix.hpp"

/*
	Alexander �ilykov, Jan 2017
*/

template <typename T>
class AbstractPreconditioner
	: virtual public AbstractMultipliableMatrix<T> {
public:
	// A := L + D + U
	// find y := [ L + w D ]^-1 . x
	virtual std::vector<T> forwSubst(std::vector<T> const & x, double w = 1.) const = 0;
	// find y := [ U + w D ]^-1 . x			
	virtual std::vector<T> backSubst(std::vector<T> const &, double w = 1.) const = 0;
	// find for y := [ D ]^-1 . x				
	virtual std::vector<T> diagSubst(std::vector<T> const &) const = 0;
	// find for y := D . x				
	virtual std::vector<T> multDiag(std::vector<T> const &) const = 0;
};