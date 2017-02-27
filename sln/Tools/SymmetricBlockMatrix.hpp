#pragma once
#include "AbstractTransposeMultipliableMatrix.hpp"

/*
	Alexander Žilyakov, Feb 2017
*/

template <typename T>
class SymmetricBlockMatrix
	: public AbstractMultipliableMatrix<T> {
	std::vector<Index> _blockStartIndicies;
	// the first element of (i, j)–block has row index _blockStartIndicies[i] and col index _blockStartIndicies[j]
	std::vector<std::vector<AbstractMultipliableMatrix<T>*>>          _diagBlocksPtr; // pointers to diag symmetric blocks
	std::vector<std::vector<AbstractTransposeMultipliableMatrix<T>*>> _lvalBlocksPtr; // pointers to blocks of lower triangular part, col by col
	T& _set(Index i, Index j) final {
		throw std::logic_error("not implemented");
	}
	T  _get(Index i, Index j) const final {
		throw std::logic_error("not implemented");
	}
public:
	SymmetricBlockMatrix(
			std::initializer_list<AbstractMultipliableMatrix<T>*>          const & diagIniList,
			std::initializer_list<AbstractTransposeMultipliableMatrix<T>*> const & lvalIniList
		) 
		: _blockStartIndicies(diagIniList.size())
		, _diagBlocksPtr(diagIniList)
		, _lvalBlocksPtr(lvalIniList) 
	{
		Index n = _diagBlocksPtr.size(), m = _lvalBlocksPtr.size();
		if ((n * n - n) / 2 != m) throw std::invalid_argument("invalid symmetric block matrix: check diag / lower triangle sizes");
		// build starting indicies
		_h = 0;
		for (Index bj = 0, k = 0; bj < n - 1; ++bj) {
			Index numbOfCols = 0; 
			if (_diagBlocksPtr[bj]) numbOfCols = _diagBlocksPtr[bj]->getOrder();
			else for (Index bi = bj + 1; bi < n; ++bi, ++k) if (_lvalBlocksPtr[k]) numbOfCols = _lvalBlocksPtr[k]->numbOfCols();
			if (_lvalBlocksPtr[k - 1]) _h = numbOfCols;
			_blockStartIndicies[bj + 1] = _blockStartIndicies[bj] + numbOfCols;
		}
		_w = (_h += _blockStartIndicies.back());
		// TODO: check block sizes
		// std::cout << "block col starting indicies: " << _blockColStartIndicies << '\n' << _blockRowStartIndicies << '\n' << _h << ' ' << _w;
	}
	SymmetricBlockMatrix& operator=(T const & val) final {
		throw std::logic_error("not implemented");
		// TODO: operator=
		//std::for_each(_blocksPtr.begin(), _blocksPtr.end(), [&](auto& row) {
		//	std::for_each(row.begin(), row.end(), [&](auto& mtx) {
		//		if (mtx) *mtx = val;
		//	});
		//});
		//return *this;
	}
	void mult(T const * by, T* result) const final {
		throw std::logic_error("not implemented");
		for (Index bi = 0; bi < _hb; ++bi) // for each block row and
			for (Index bj = 0; bj < _wb; ++bj) // col, do…
				if (_blocksPtr[bi][bj]) _blocksPtr[bi][bj]->mult(by + _blockColStartIndicies[bj], result + _blockRowStartIndicies[bi]);
	}
};
