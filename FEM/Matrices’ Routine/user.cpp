#include <fstream>
#include "CRSMatrix.hpp"

int main() {
	std::ifstream inputMatrix("matrix.txt");
	size_t order, nonzeros;
	inputMatrix >> order >> nonzeros;
	CRSMatrix A(order, nonzeros);
	inputMatrix >> A;
	// mult by vector
	std::vector<REAL> iVec(order, 1.);
	std::cout << A * iVec;
	// modify and print
	A(1, 1) = 56.;
	A.print(std::cout); // trad form
	std::cout << A; // sparse form

	// TODO:
	// std::vector<REAL> x(order), b(order, 0.);
	// b[0] = 1.;
	// b[order - 1] = 1.;
	// x = A.CG(b); // or x = b / A
	
	return 0;
}