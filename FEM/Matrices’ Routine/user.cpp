#include <fstream>
#include "CRSMatrix.hpp"
#include "BandMatrix.hpp"
#include "SymmetricMatrix.hpp"

int main() {
	SymmetricMatrix I(3);
	try {
		std::cout << I << '\n';
	}
	catch (std::exception& e) {
		std::cout << e.what();
	}

	std::ifstream inputMatrixCRS("matrix_crs.txt"),
				  inputMatrixBand("matrix_band.txt");
	size_t order, nonzeros, bandWidth;
	// CRS
	inputMatrixCRS >> order >> nonzeros;
	CRSMatrix A(order, nonzeros);
	inputMatrixCRS >> A;
	// mult by vector
	std::vector<REAL> iVec(order, 1.);
	std::cout << A * iVec << '\n';
	// modify and print
	A(1, 1) = 56.;
	A.print(); // trad form
	std::cout << A << '\n'; // sparse form
	// Band
	inputMatrixBand >> order >> bandWidth;
	BandMatrix B(order, bandWidth);
	inputMatrixBand >> B;
	std::vector<REAL> x(order), f(order, 0.);
	f[0] = 1.;
	f[order - 1] = 1.;
	std::cout << f / B; // or x = B.solve(f)
	return 0;
}