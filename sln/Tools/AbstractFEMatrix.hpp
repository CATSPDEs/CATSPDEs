#pragma once
#include "AbstractSparseMatrix.hpp"
#include "AdjacencyList.hpp"

/*
	Alexander Žilyakov, Sep 2016
*/

template <typename T>
class AbstractFEMatrix
	: virtual public AbstractSparseMatrix<T> {
public:
	// construct matrix pattern from adjacency list 
	// of finite element mesh verticies
	virtual AbstractFEMatrix& generatePatternFrom(AdjacencyList const &) = 0;
};