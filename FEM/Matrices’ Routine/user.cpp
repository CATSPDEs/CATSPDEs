#include <fstream>
#include "vector.hpp" // vector operations (*, +, ...)
#include "CRSMatrix.hpp"
#include "BandMatrix.hpp"

int main() {
	std::ifstream inputMatrixCRS("matrix_crs.txt"),
				  inputMatrixBand("matrix_band.txt");
	size_t order, nonzeros, bandWidth;
	// CRS
	inputMatrixCRS >> order >> nonzeros;
	CRSMatrix A(order, nonzeros);
	inputMatrixCRS >> A;
	// mult by vector
	std::vector<double> iVec(order, 1.);
	std::cout << A * iVec << '\n';
	// modify and print
	A(1, 1) = 56.;
	A.save(std::cout); // trad form
	std::cout << A << '\n'; // sparse form
	// Band
	inputMatrixBand >> order >> bandWidth;
	BandMatrix B(order, bandWidth);
	inputMatrixBand >> B;
	std::vector<double> x(order), f(order, 0.);
	f[0] = 1.;
	f[order - 1] = 1.;
	std::cout << f / B << std::endl; // or x = B.solve(f)
	B(1, 0) = 56.;
	B.save(std::cout); // (almost) LU-decomposition is now stored in B 
	return 0;
}