#include <fstream>
#include <string>
#include "CSRMatrix.hpp"
#include "CSCMatrix.hpp"

int main() {
	SingletonLogger& logger = SingletonLogger::instance();
	logger.beg("choose matrix format from CATSPDEs collection");
		string matrixType;
		vector<string> availableTypes = { "CSR", "CSC" };
		logger.log("available matrix types");
		for (auto const & t : availableTypes)
			logger.mes("* " + t);
		while (find(availableTypes.begin(), availableTypes.end(), matrixType) == availableTypes.end()) {
			logger.inp("enter matrix type");
			cin >> matrixType;
		}
	logger.end();
	try {
		/*
			CRS–matrix
		*/
		if (matrixType == "CSR") {
			logger.beg("CSR matrix test");
				ifstream iCSR("Mathematica/matrices/csr.dat"),
						 iU("Mathematica/matrices/csr_vector.dat"),
						 iV("Mathematica/matrices/csr_vector_t.dat");
				ofstream oCSRxU("Mathematica/matrices/csr_mult.dat"),
						 oCSRxV("Mathematica/matrices/csr_mult_t.dat");
				size_t rows, cols, nnz;
				iCSR >> rows >> cols >> nnz;
				CSRMatrix<double> CSR(rows, cols, nnz);
				iCSR >> CSR;
				vector<double> u(cols), v(rows);
				iU >> u;
				iV >> v;
				oCSRxU << CSR * u;
				oCSRxV << CSR.t() * v;
			logger.end();
		}
		/*
			CSC–matrix
		*/
		else if (matrixType == "CSC") {
			logger.beg("CSC matrix test");
				logger.beg("mult() / multTranspose() test");
					ifstream iCSC("Mathematica/matrices/csc.dat"),
							 iU("Mathematica/matrices/csc_vector.dat"),
							 iV("Mathematica/matrices/csc_vector_t.dat");
					ofstream oCSCxU("Mathematica/matrices/csc_mult.dat"),
							 oCSCxV("Mathematica/matrices/csc_mult_t.dat");
					size_t rows, cols, nnz;
					iCSC >> rows >> cols >> nnz;
					CSCMatrix<double> CSC(rows, cols, nnz);
					iCSC >> CSC;
					vector<double> u(cols), v(rows);
					iU >> u;
					iV >> v;
					oCSCxU << CSC * u;
					oCSCxV << CSC.t() * v;
				logger.end();
				logger.beg("Harwell-Boeing i/o test");
					HBMatrix<double>* HB = nullptr;
					HBMatrix<double>::loadHarwellBoeing("Mathematica/HarwellBoeing/input/beacxc.rra", HB, logger);
					//readHarwellBoeingStruct();
					//HB.save(oHBDense);
				logger.end();
			logger.end();
		}
	}
	catch (exception const & e) {
		logger.err(e.what());
	}
	return 0;
}