#pragma once
#include "AbstractMultipliableMatrix.hpp"

/*
	Alexander Žilyakov, Feb 2017
*/

template <typename T>
class BlockMatrix
	: public AbstractMultipliableMatrix<T> { // aggregation
	std::vector<Index> _blockRowStartIndicies, _blockColStartIndicies;
	// the first element of (i, j)–block has row index _blockRowStartIndicies[i] and col index _blockColStartIndicies[j]
	std::vector<std::vector<AbstractMultipliableMatrix<T>*>> _blocksPtr; // pointers to blocks
	T& _set(Index i, Index j) final {
		throw std::logic_error("not implemented");
	}
	T  _get(Index i, Index j) const final {
		throw std::logic_error("not implemented");
	}
public:
	BlockMatrix(
		std::initializer_list<std::initializer_list<AbstractMultipliableMatrix<T>*>> const & iniList
	) 
		: _blockRowStartIndicies(iniList.size())
		, _blockColStartIndicies(iniList.begin()->size()) 
	{
		Index n = _blockRowStartIndicies.size(), m = _blockColStartIndicies.size();
		_blocksPtr.reserve(n * m);
		for (auto const & row : iniList) {
			if (row.size() != m) throw std::invalid_argument("ini list does not represent a block matrix");
			_blocksPtr.emplace_back(row);
		}
		std::vector<Index> blockRowSizes(n), blockColSizes(m);
		auto modifyOnce = [](Index& what, Index val) {
			if (what && what != val) throw std::invalid_argument("invalid block sizes of block matrix");
			what = val;
		};
		for (Index bi = 0; bi < n; ++bi)
			for (Index bj = 0; bj < m; ++bj)
				if (_blocksPtr[bi][bj]) {
					modifyOnce(blockRowSizes[bi], _blocksPtr[bi][bj]->numbOfRows());
					modifyOnce(blockColSizes[bj], _blocksPtr[bi][bj]->numbOfCols());
				}
		for (Index bi = 1; bi < n; ++bi) _blockRowStartIndicies[bi] = _blockRowStartIndicies[bi - 1] + blockRowSizes[bi - 1];
		_h = _blockRowStartIndicies.back() + blockRowSizes.back();
		for (Index bj = 1; bj < m; ++bj) _blockColStartIndicies[bj] = _blockColStartIndicies[bj - 1] + blockColSizes[bj - 1];
		_w = _blockColStartIndicies.back() + blockColSizes.back();
		// std::cout << "block col / row starting indicies: " << _blockColStartIndicies << '\n' << _blockRowStartIndicies << '\n' << _h << ' ' << _w;
	}
	BlockMatrix& operator=(T const & val) final {
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
		Index n = _blockRowStartIndicies.size(), m = _blockColStartIndicies.size();
		for (Index bi = 0; bi < n; ++bi) // for each block row and
			for (Index bj = 0; bj < m; ++bj) // col, do…
				if (_blocksPtr[bi][bj]) _blocksPtr[bi][bj]->mult(by + _blockColStartIndicies[bj], result + _blockRowStartIndicies[bi]);
	}
};