#include <iostream>
#include "CRSMatrix.h"

int main() {
	std::ifstream inputMatrix("matrix.txt");
	CRSMatrix A(inputMatrix);
	unsigned i, n = A.getOrder();
	std::vector<REAL> u(n), v(n);
	for (i = 0; i < n; ++i)
		u[i] = i + 1;
	v = A * u;
	A.print();
	return 0;
}