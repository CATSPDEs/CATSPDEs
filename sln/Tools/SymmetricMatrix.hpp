#pragma once
#include "AbstractSparseMatrix.hpp"

template <typename T>
class SymmetricMatrix : public AbstractSparseMatrix<T> {
	std::vector<T> _val; // := upper triangle, row by row
	// virtual methods to be implemented
	T& _set(Index i, Index j) final {
		if (i > j) std::swap(i, j);
		return _val[j + _w * i - (i * i + i) / 2];
	}
	T _get(Index i, Index j) const final {
		if (i > j) std::swap(i, j);
		return _val[j + _w * i - (i * i + i) / 2];
	}
public:
	explicit SymmetricMatrix(Index n = 1) : AbstractMatrix<T>(n, n), _val(nnz()) {}
	SymmetricMatrix(std::initializer_list<T> const & iniList) // create matrix from 
		: AbstractMatrix<T>(round(.5 * (sqrt(1. + 8. * iniList.size()) - 1.))) 
		, _val(iniList)
		{}
	// virtual methods to be implemented
	Index nnz() const final { return (_w * _w + _w) / 2; }
	SymmetricMatrix& operator=(T const & val) final {
		std::fill(_val.begin(), _val.end(), val);
		return *this;
	}
	SymmetricMatrix& importSparse(std::istream& from = cin) final {
		from >> _val;
		return *this;
	}
	void exportSparse(std::ostream& to = cout) const final {
		to << _w << '\n' << _val;
	}
	// almost same methods from base class (in order to work w/ strings instead of streams)
	using AbstractSparseMatrix::importSparse;
	using AbstractSparseMatrix::exportSparse;
	SymmetricMatrix& operator*=(T const & val) final {
		_val *= val;
		return *this;
	}
};