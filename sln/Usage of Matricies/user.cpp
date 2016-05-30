#include <fstream>
#include "CSRMatrix.hpp"
#include "BandMatrix.hpp"
#include "DenseMatrix.hpp"
#include "SymmetricCSlRMatrix.hpp"
#include "krylov.hpp"

int main() {
	std::ifstream inA("testCG_A.txt");
	std::ifstream inb("testCG_b.txt");
	SymmetricCSlRMatrix A_cslr(4, 3);
	A_cslr.loadSparse(inA);
	A_cslr.save(std::cout);
	std::vector<double> b_cslr(4);
	inb >> b_cslr;
	std::vector<double> x0(4, 0);
	std::cout << CG(A_cslr, b_cslr, x0, 10e-7);
	std::cout << std::endl << std::endl;
	std::ifstream inputMatrixCRS("matrix_crs.txt"),
				  inputMatrixBand("matrix_band.txt");
	size_t order, nonzeros, bandWidth;
	// CRS
	inputMatrixCRS >> order >> nonzeros;
	CSRMatrix A(order, nonzeros);
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