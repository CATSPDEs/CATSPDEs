#include "AbstractSquareMatrix.hpp"

AbstractSquareMatrix::AbstractSquareMatrix(size_t n) : _n(n) {
	if (_n < 1) throw std::out_of_range("order of matrix must be at least one");
}

void AbstractSquareMatrix::print() {
	size_t i, j;
	if (sizeof(REAL) == 4) std::cout.precision(6);
	else std::cout.precision(14);
	std::cout << std::scientific;
	for (i = 0; i < _n; ++i) {
		for (j = 0; j < _n; ++j) 
			std::cout << get(i, j) << ' ';
		std::cout << '\n';
	}
}

std::istream& operator>>(std::istream& from, AbstractSquareMatrix& A) { return A.load(from); }
std::ostream& operator<<(std::ostream& to, AbstractSquareMatrix const & A) { return A.save(to); }
std::vector<REAL> operator*(AbstractSquareMatrix const & A, std::vector<REAL> const & b) { return A.mult(b); }
std::vector<REAL> operator/(std::vector<REAL> const & b, AbstractSquareMatrix & A) { return A.solve(b); }