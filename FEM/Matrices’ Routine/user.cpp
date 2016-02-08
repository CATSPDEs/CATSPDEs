#include "SLAE.h"

int main() {
	std::ifstream inputMatrix("matrix.txt");
	CRSMatrix A(inputMatrix);
	unsigned i, n = A.getOrder();
	A(1, 1) = 56.;
	//std::vector<REAL> x(n), b(n, 0.);
	//b[0] = 1.;
	//b[n - 1] = 1.;
	//SLAE system(A, b);
	//x = system.solve();
	A.print();
	return 0;
}