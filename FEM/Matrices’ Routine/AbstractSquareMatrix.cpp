#include "AbstractSquareMatrix.hpp"

std::ostream& AbstractSquareMatrix::print(std::ostream& output) const {
	size_t i, j;
	if (sizeof(REAL) == 4) output.precision(6);
	else output.precision(14);
	output << std::scientific;
	for (i = 0; i < _n; ++i) {
		for (j = 0; j < _n; ++j) output << (*this)(i, j) << ' ';
		output << '\n';
	}
	output << std::endl;
	return output;
}

std::istream& operator>>(std::istream& from, AbstractSquareMatrix& A) { return A.load(from); }
std::ostream& operator<<(std::ostream& to, AbstractSquareMatrix const & A) { return A.save(to); }
std::vector<REAL> operator*(AbstractSquareMatrix const & A, std::vector<REAL> const & b) { return A.mult(b); }
std::vector<REAL> operator/(std::vector<REAL> const & b, AbstractSquareMatrix const & A) { return A.solve(b); }