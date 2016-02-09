#include <fstream>
#include "CRSMatrix.h"

int main() {
	std::ifstream inputMatrix("matrix.txt");
	size_t order, nonzeros;
	inputMatrix >> order >> nonzeros;
	CRSMatrix A(order, nonzeros);
	inputMatrix >> A;
	A(1, 1) = 56.;
	std::cout << A;
	// TODO:
	// std::vector<REAL> x(order), b(order, 0.);
	// b[0] = 1.;
	// b[order - 1] = 1.;
	// x = A.CG(b); // or x = b / A
	return 0;
}