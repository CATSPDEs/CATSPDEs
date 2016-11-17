#pragma once
#include "AbstractSparseMatrix.hpp"

template <typename T>
class SymmetricContainer : public AbstractSparseMatrix<T> {
	std::vector<T> _val;
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
	explicit SymmetricContainer(Index n = 1) : AbstractMatrix<T>(n, n), _val(nnz()) {}
	~SymmetricContainer() {}
	// virtual methods to be implemented
	Index nnz() const final { return (_w * _w + _w) / 2; }
	SymmetricContainer& operator=(T const & val) final {
		std::fill(_val.begin(), _val.end(), val);
		return *this;
	}
	SymmetricContainer& importSparse(std::istream& from = cin) final {
		from >> _val;
		return *this;
	}
	void exportSparse(std::ostream& to = cout) const final {
		to << _w << '\n' << _val;
	}
};