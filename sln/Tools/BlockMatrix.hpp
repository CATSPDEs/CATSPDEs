#pragma once
#include "AbstractMultipliableMatrix.hpp"

/*
	Alexander Žilyakov, Feb 2017
*/

template <typename T>
class BlockMatrix
	: public AbstractMultipliableMatrix<T> {
	Index _hb, _wb; // numb of row and col blocks
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
	BlockMatrix(std::initializer_list<std::initializer_list<AbstractMultipliableMatrix<T>*>> const & iniList) 
		: _hb(iniList.size())
		, _wb(iniList.begin()->size())
		, _blockRowStartIndicies(_hb) 
		, _blockColStartIndicies(_wb) {
		for (auto const & row : iniList) {
			if (row.size() != _wb) throw std::invalid_argument("ini list does not represent a block matrix");
			_blocksPtr.emplace_back(row);
		}
		Index numbOfCols, numbOfRows;
		// (1) build blocks’ col starting indicies
		for (Index bj = 0; bj < _wb - 1; ++bj) {
			_blockColStartIndicies[bj + 1] = _blockColStartIndicies[bj];
			numbOfCols = 0;
			for (Index bi = 0; bi < _hb; ++bi) 
				if (_blocksPtr[bi][bj]) {
					if (!numbOfCols) numbOfCols = _blocksPtr[bi][bj]->numbOfCols();
					else if (numbOfCols != _blocksPtr[bi][bj]->numbOfCols()) throw std::invalid_argument("blocks of different width: check block col #" + std::to_string(bj));
				}
			_blockColStartIndicies[bj + 1] += numbOfCols;
		}
		// (2) check numb of cols of last block col for equality and compute _w
		numbOfCols = 0;
		for (Index bi = 0; bi < _hb; ++bi)
			if (_blocksPtr[bi][_wb - 1]) {
				if (!numbOfCols) numbOfCols = _blocksPtr[bi][_wb - 1]->numbOfCols();
				else if (numbOfCols != _blocksPtr[bi][_wb - 1]->numbOfCols()) throw std::invalid_argument("blocks of different width: check last block col");
			}
		_w = _blockColStartIndicies.back() + numbOfCols;
		// ″ for row indicies
		for (Index bi = 0; bi < _hb - 1; ++bi) {
			_blockRowStartIndicies[bi + 1] = _blockRowStartIndicies[bi];
			numbOfRows = 0;
			for (Index bj = 0; bj < _wb; ++bj)
				if (_blocksPtr[bi][bj]) {
					if (!numbOfRows) numbOfRows = _blocksPtr[bi][bj]->numbOfRows();
					else if (numbOfRows != _blocksPtr[bi][bj]->numbOfRows()) throw std::invalid_argument("blocks of different height: check block row #" + std::to_string(bi));
				}
			_blockRowStartIndicies[bi + 1] += numbOfRows;
		}
		numbOfRows = 0;
		for (Index bj = 0; bj < _wb; ++bj)
			if (_blocksPtr[_hb - 1][bj]) {
				if (!numbOfRows) numbOfRows = _blocksPtr[_hb - 1][bj]->numbOfRows();
				else if (numbOfRows != _blocksPtr[_hb - 1][bj]->numbOfRows()) throw std::invalid_argument("blocks of different height: check last block row");
			}
		_h = _blockRowStartIndicies.back() + numbOfRows;
		// std::cout << "block col starting indicies: " << _blockColStartIndicies << '\n' << _blockRowStartIndicies << '\n' << _h << ' ' << _w;
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
		for (Index bi = 0; bi < _hb; ++bi) // for each block row and
			for (Index bj = 0; bj < _wb; ++bj) // col, do…
				if (_blocksPtr[bi][bj]) _blocksPtr[bi][bj]->mult(by + _blockColStartIndicies[bj], result + _blockRowStartIndicies[bi]);
	}
};

//template <typename T>
//class SaddlePointMatrix
//	: public AbstractMultipliableMatrix<T> {
//	T& _set(Index i, Index j) final {
//		auto n = A11.getOrder();
//		// A11 
//		if (i < n && j < n) return A11(i, j);
//		// A22
//		if (n <= i && i < 2 * n &&
//			n <= j && j < 2 * n) return A11(i - n, j - n);
//		// to lower triangle
//		if (i < j) std::swap(i, j);
//		// not–A
//		if (2 * n <= i) {
//			// B11
//			if (j < n) return B1(i - 2 * n, j);
//			// B22
//			if (j <= 2 * n) return B2(i - 2 * n, j - 2 * n);
//		}
//		throw std::invalid_argument("pattern of sparse matrix does not allow to change element w/ these indicies");
//	}
//	T  _get(Index i, Index j) const final {
//		auto n = A11.getOrder();
//		// A11 
//		if (i < n && j < n) return A11(i, j);
//		// A22
//		if (n <= i && i < 2 * n &&
//			n <= j && j < 2 * n) return A11(i - n, j - n);
//		// to lower triangle
//		if (i < j) std::swap(i, j);
//		// not–A
//		if (2 * n <= i) {
//			// B11
//			if (j < n) return B1(i - 2 * n, j);
//			// B22
//			if (j < 2 * n) return B2(i - 2 * n, j - n);
//		}
//		return 0.;
//	}
//public:
//	//
//	//  | A   B^T |
//	//  | B   0   |
//	//
//	// A:
//	//
//	//  | A11   0   |
//	//  | 0     A22 |
//	//
//	// A11 = A22
//	//
//	// B:
//	//
//	//  | B1  B2 |
//	//
//
//	CSlCMatrix<T> A11;
//	CSCMatrix<double> B1, B2;
//	SaddlePointMatrix(CSlCMatrix<T> const & a11, CSCMatrix<double> const & b1, CSCMatrix<double> const & b2) 
//		: A11(a11), B1(b1), B2(b2), AbstractMatrix<T>(2 * a11.getOrder() + b1.numbOfRows())
//	{}
//	SaddlePointMatrix& operator=(T const & val) final {
//		A11 = val; B1 = val; B2 = val;
//		return *this;
//	}
//	void mult(T const * by, T* result) const final {
//		// ...
//	}
//};