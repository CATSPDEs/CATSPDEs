#include <fstream>
#include <string>
#include "CSRMatrix.hpp"
#include "CSCMatrix.hpp"

int main() {
	string matrixType;
	vector<string> availableTypes = { "CSR", "CSC" };
	cout << "available matrix types:\n";
	for (auto const & t : availableTypes)
		cout << "    * " << t << '\n';
	while (find(availableTypes.begin(), availableTypes.end(), matrixType) == availableTypes.end()) {
		cout << "enter matrix type: ";
		cin >> matrixType;
	}
	try {
		if (matrixType == "CSR") {
			ifstream iCSR("Mathematica/matrices/csr.dat"),
				iU("Mathematica/matrices/csr_vector.dat"),
				iV("Mathematica/matrices/csr_vector_t.dat");
			ofstream oCSRxU("Mathematica/matrices/csr_mult.dat"),
				oCSRxV("Mathematica/matrices/csr_mult_t.dat");
			size_t rows, cols, nnz;
			iCSR >> rows >> cols >> nnz;
			CSRMatrix<double> CSR(rows, cols, nnz);
			iCSR >> CSR;
			CSR.saveSparse();
			vector<double> u(cols), v(rows);
			iU >> u;
			iV >> v;
			oCSRxU << CSR * u;
			oCSRxV << CSR.t() * v;
		}
		else if (matrixType == "CSC") {
			ifstream iCSC("Mathematica/matrices/csc.dat"),
			         iU("Mathematica/matrices/csc_vector.dat"),
			         iV("Mathematica/matrices/csc_vector_t.dat");
			ofstream oCSCxU("Mathematica/matrices/csc_mult.dat"),
			         oCSCxV("Mathematica/matrices/csc_mult_t.dat");
			size_t rows, cols, nnz;
			iCSC >> rows >> cols >> nnz;
			CSCMatrix<double> CSC(rows, cols, nnz);
			iCSC >> CSC;
			CSC.save();
			vector<double> u(cols), v(rows);
			iU >> u;
			iV >> v;
			oCSCxU << CSC * u;
			oCSCxV << CSC.t() * v;
			// HB
			ifstream iHB("HarwellBoeing/beacxc.rra");
			ofstream oHBDense("Mathematica/beacxc.dat");
			HBMatrix<double> HB(iHB);
			HB.save(oHBDense);
		}
	}
	catch (exception const & e) {
		cout << "error: " << e.what();
	}
	return 0;
}