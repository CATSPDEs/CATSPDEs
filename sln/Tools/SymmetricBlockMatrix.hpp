#pragma once
#include "AbstractTransposeMultipliableMatrix.hpp"

/*
	Alexander Žilyakov, Feb 2017
*/

enum class BlockSymmetryType { symmetric, antisymmetric };

template <typename T>
class SymmetricBlockMatrix 
	: public AbstractMultipliableMatrix<T> { // aggregation
	short _bsign; // + for symmetric, – for antisymmetric
	std::vector<Index> _blockStartIndicies;
	// the first element of (i, j)–block has row index _blockStartIndicies[i] and col index _blockStartIndicies[j]
	std::vector<AbstractMultipliableMatrix<T>*>          _diagBlocksPtr; // pointers to diag symmetric blocks
	std::vector<AbstractTransposeMultipliableMatrix<T>*> _lvalBlocksPtr; // pointers to blocks of lower triangular part, col by col
	T& _set(Index i, Index j) final {
		throw std::logic_error("not implemented");
	}
	T  _get(Index i, Index j) const final {
		throw std::logic_error("not implemented");
	}
public:
	SymmetricBlockMatrix(
			std::initializer_list<AbstractMultipliableMatrix<T>*>          const & diagIniList,
			std::initializer_list<AbstractTransposeMultipliableMatrix<T>*> const & lvalIniList,
			BlockSymmetryType stype = BlockSymmetryType::symmetric
		) 
		: _bsign(stype == BlockSymmetryType::symmetric ? 1 : -1)
		, _blockStartIndicies(diagIniList.size())
		, _diagBlocksPtr(diagIniList)
		, _lvalBlocksPtr(lvalIniList) 
	{
		Index n = _diagBlocksPtr.size(), m = _lvalBlocksPtr.size();
		if ((n * n - n) / 2 != m) throw std::invalid_argument("invalid symmetric block matrix: check numb of blocks in diag or lower triangular part");
		std::vector<Index> blockSizes(n);
		auto modifyOnce = [](Index& what, Index val) {
			if (what && what != val) throw std::invalid_argument("invalid block sizes of symmetric block matrix");
			what = val;
		};
		for (Index bj = 0, k = 0; bj < n; ++bj) {
			if (_diagBlocksPtr[bj]) modifyOnce(blockSizes[bj], _diagBlocksPtr[bj]->getOrder());
			for (Index bi = bj + 1; bi < n; ++bi, ++k)
				if (_lvalBlocksPtr[k]) modifyOnce(blockSizes[bi], _lvalBlocksPtr[k]->numbOfRows());
		}
		for (Index bj = 1; bj < n; ++bj)
			_blockStartIndicies[bj] = _blockStartIndicies[bj - 1] + blockSizes[bj - 1];
		_h = _w = _blockStartIndicies.back() + blockSizes.back();
		// std::cout << "block starting indicies: " << _blockStartIndicies << '\n' << '\n' << _h;
	}
	SymmetricBlockMatrix& operator=(T const & val) final {
		throw std::logic_error("not implemented");
	}
	void mult(T const * by, T* result) const final {
		Index n = _diagBlocksPtr.size();
		for (Index bj = 0, k = 0; bj < n; ++bj) {
			if (_diagBlocksPtr[bj]) 
				_diagBlocksPtr[bj]->mult(by + _blockStartIndicies[bj], result + _blockStartIndicies[bj]);
			for (Index bi = bj + 1; bi < n; ++bi, ++k)
				if (_lvalBlocksPtr[k]) {
					_lvalBlocksPtr[k]->mult           (by + _blockStartIndicies[bj], result + _blockStartIndicies[bi]);

					_lvalBlocksPtr[k]->multByTranspose(by + _blockStartIndicies[bi], result + _blockStartIndicies[bj]);
				}
		}
	}
};
