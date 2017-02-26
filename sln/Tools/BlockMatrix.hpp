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
	// the first element of (i, j)–block has row index _blockRowStartIndicies[i] and _blockColStartIndicies[j]
	std::vector<std::vector<AbstractMultipliableMatrix<T>*>> _A; // pointers to blocks

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
			_A.emplace_back(row);
		}
		// (1) build blocks’ col starting indicies
		for (Index bj = 0; bj < _wb - 1; ++bj) {
			_blockColStartIndicies[bj + 1] = _blockColStartIndicies[bj];
			Index numbOfCols = 0;
			for (Index bi = 0; bi < _hb; ++bi) 
				if (_A[bi][bj]) {
					if (!numbOfCols) numbOfCols = _A[bi][bj]->numbOfCols();
					else if (numbOfCols != _A[bi][bj]->numbOfCols()) throw std::invalid_argument("blocks of different width: check block col #" + std::to_string(bj));
				}
			_blockColStartIndicies[bj + 1] += numbOfCols;
		}
		// (2) check numb of cols of last block col for equality 
		for (Index bi = 0, numbOfCols = 0; bi < _hb; ++bi)
			if (_A[bi][_wb - 1]) {
				if (!numbOfCols) numbOfCols = _A[bi][_wb - 1]->numbOfCols();
				else if (numbOfCols != _A[bi][_wb - 1]->numbOfCols()) throw std::invalid_argument("blocks of different width: check last block col");
			}
		// ″ for row indicies
		for (Index bi = 0; bi < _hb - 1; ++bi) {
			_blockRowStartIndicies[bi + 1] = _blockRowStartIndicies[bi];
			Index numbOfRows = 0;
			for (Index bj = 0; bj < _wb; ++bj)
				if (_A[bi][bj]) {
					if (!numbOfRows) numbOfRows = _A[bi][bj]->numbOfRows();
					else if (numbOfRows != _A[bi][bj]->numbOfRows()) throw std::invalid_argument("blocks of different height: check block row #" + std::to_string(bi));
				}
			_blockRowStartIndicies[bi + 1] += numbOfRows;
		}
		for (Index bj = 0, numbOfRows = 0; bj < _wb; ++bj)
			if (_A[_hb - 1][bj]) {
				if (!numbOfRows) numbOfRows = _A[_hb - 1][bj]->numbOfRows();
				else if (numbOfRows != _A[_hb - 1][bj]->numbOfRows()) throw std::invalid_argument("blocks of different height: check last block row");
			}
		std::cout << "block col starting indicies: " << _blockColStartIndicies << '\n' << _blockRowStartIndicies;
	}
	BlockMatrix& operator=(T const & val) final {
		throw std::logic_error("not implemented");
		// TODO: operator=
		//std::for_each(_A.begin(), _A.end(), [&](auto& row) {
		//	std::for_each(row.begin(), row.end(), [&](auto& mtx) {
		//		if (mtx) *mtx = val;
		//	});
		//});
		//return *this;
	}
	void mult(T const * by, T* result) const final {
		throw std::logic_error("not implemented");
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