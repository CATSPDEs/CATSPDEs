#include <fstream>
#include "CRSMatrix.hpp"
#include "BandMatrix.hpp"
#include "DenseMatrix.hpp"

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
	try {
		A(3, 1) = 56.;
	}
	catch (std::exception const & e) {
		std::cout << e.what();
	}
	A.save(); // trad form
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
	B.save(); // (almost) LU-decomposition is now stored in B 
	// Dense
	DenseMatrix D(2);
	D(0, 0) = 1;
	D(1, 1) = 2;
	D(0, 1) = D(1, 0) = 0;
	D.save();
	std::vector<double> b = {1, 2};
	std::cout << b / D << std::endl;
	return 0;
}