#pragma once
#include"AbstractSparseMatrix.hpp"
#include"IRealMatrix.hpp"

class SkylineMatrix : public AbstractSparseMatrix<double>, public IRealMatrix {
	std::vector<size_t> _ptr,
						_col;
	std::vector<double> _val,
						_diag;
	size_t _nnz() const { return _ptr[_n] + _n; }
	double& _set(size_t, size_t);
	double _get(size_t, size_t) const;
public:
	SkylineMatrix(size_t, size_t);
	~SkylineMatrix() {}
	std::vector<double> solve(std::vector<double> const &);
	std::vector<double> mult(std::vector<double> const &) const;
	std::istream& loadSparse(std::istream&);
	std::ostream& saveSparse(std::ostream&) const;
	friend SkylineMatrix operator*(SkylineMatrix const &, SkylineMatrix const &);
};