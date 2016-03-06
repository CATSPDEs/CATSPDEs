#include "AbstractRealMatrix.hpp"

void AbstractRealMatrix::print() {
	size_t i, j;
	std::cout.precision(14);
	std::cout << std::scientific;
	for (i = 0; i < _n; ++i) {
		for (j = 0; j < _n; ++j)
			std::cout << get(i, j) << ' ';
		std::cout << '\n';
	}
}

std::istream& operator>>(std::istream& from, AbstractRealMatrix& A) { return A.load(from); }
std::ostream& operator<<(std::ostream& to, AbstractRealMatrix const & A) { return A.save(to); }
std::vector<double> operator*(AbstractRealMatrix const & A, std::vector<double> const & b) { return A.mult(b); }
std::vector<double> operator/(std::vector<double> const & b, AbstractRealMatrix & A) { return A.solve(b); }