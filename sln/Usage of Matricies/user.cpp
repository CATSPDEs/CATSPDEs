#include <complex> // for complex–valued matrices
#include "CSCMatrix.hpp"
#include "SymmetricCSlCMatrix.hpp"
#include "SymmetricMatrix.hpp"
#include "DenseMatrix.hpp"
#include "SingletonLogger.hpp"

/*
	Alexander Žilykov, Aug 2016
*/

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::complex;

int main() {
	// common vars
	string iPath, oPath;
	vector<string> HarwellBoeingMatricies = {
		"mathematica.rra (real    rectangular assembled) -- simple artificial example",
		"mathematica.rsa (real    symmetric   assembled) -- simple artificial example",
		"illc1033.rra    (real    rectangular assembled)",
		"e40r5000.rua    (real    unsymmetric assembled)",
		"qc2534.cua      (complex unsymmetric assembled)",
		"young1c.csa     (complex symmetric   assembled)",
		"cegb2802.pse    (pattern symmetric   elemental)"
	};
	// logger
	SingletonLogger& logger = SingletonLogger::instance();
	try {
		auto testNum = logger.opt("choose matrix format from CATSPDEs collection", { 
			"CSC", 
			"CSlC", 
			"SymmetricCSlC", 
			"SymmetricMatrix",
			"DenseMatrix"
		});
		if (testNum == 0) {
			logger.beg("CSC matrix test");
				logger.beg("mult() / multTranspose() test");
					// i/o file pathes
					iPath = "Mathematica/output/CSC/";
					oPath = "output/CSC/";
					logger.beg("read sparse matrix A, vectors u and v from " + iPath);
						CSCMatrix<double> A; 
						A.importSparse(iPath + "A.dat"); // rect CSC matrix in CATSPDEs sparse format
						vector<double> u(A.numbOfCols()), v(A.numbOfRows());
						import(u, iPath + "u.dat"); // vector for mult() test
						import(v, iPath + "v.dat"); // vector for multTranspose() test
					logger.end();
					if (logger.yes("print A in dense, u, and v")) {
						logger.buf << "A in dense format:\n";
						A.export(logger.buf);
						logger.buf << "\nu = " << u << "\nv = " << v;
						logger.log();
					}
					logger.beg("compute multiplications A * u, A^T * v\nsave results to " + oPath);
						auto AxU = A * u;
						auto AxV = A.t() * v;
						export(AxU, oPath + "AxU.dat");
						export(AxV, oPath + "AxV.dat");
					logger.end();
					if (logger.yes("print results of multiplications")) {
						logger.buf << "A * u   = " << AxU << '\n' 
						           << "A^T * v = " << AxV;
						logger.log();
					}
					logger.log("check results w/ Mathematica/CSC.nb!");
				logger.end();
				logger.beg("Harwell-Boeing i/o test");
					iPath = "HarwellBoeing/";
					oPath = "output/CSC/HarwellBoeing/";
					auto iMatrixNum = logger.opt("choose input matrix", HarwellBoeingMatricies);
					// leave only *.rra, *.rua etc.
					for_each(HarwellBoeingMatricies.begin(), HarwellBoeingMatricies.end(), [](string& s) { s.resize(s.find_first_of(' ')); });
					string iHBPath = iPath + HarwellBoeingMatricies[iMatrixNum],
					       oHBPath = oPath + HarwellBoeingMatricies[iMatrixNum];
					if (logger.opt("choose T for CSCMatrix<T> B", { "double", "complex<double>" }) == 0) { // real
						logger.beg("load info from " + iHBPath);		
							CSCMatrix<double> B;
							auto header = B.importHarwellBoeing(iHBPath);
							logger.buf << "HB header:\n" << header;
							logger.log();
						logger.end();
						logger.beg("save matrix to " + oHBPath);
							B.exportHarwellBoeing(oHBPath, { { "title", "My real HB matrix" } });
						logger.end();
					}
					else { // complex
						logger.beg("load info from " + iHBPath);
							CSCMatrix<complex<double>> B;
							auto header = B.importHarwellBoeing(iHBPath);
							logger.buf << "HB header:\n" << header;
							logger.log();
						logger.end();
						logger.beg("save matrix to " + oHBPath);
							B.exportHarwellBoeing(oHBPath, { { "title", "My real HB matrix" } });
						logger.end();
					}
				logger.end();
			logger.end();
		}
		else if (testNum == 2) {
			logger.beg("Symmetric CSlC matrix test");
				logger.beg("mult() test");
					// i/o file pathes
					iPath = "Mathematica/output/SymmetricCSlC/";
					oPath = "output/SymmetricCSlC/";
					logger.beg("read sparse matrix A and vector u from " + iPath);
						SymmetricCSlCMatrix<double> A;
						A.importSparse(iPath + "A.dat");
						vector<double> u(A.getOrder());
						import(u, iPath + "u.dat");
					logger.end();
					if (logger.yes("print A (in dense format) and u")) {
						logger.buf << "A in dense format:\n";
						A.export(logger.buf);
						logger.buf << "u = " << u;
						logger.log();
					}
					logger.beg("compute A * u and save results to " + oPath);
						auto AxU = A * u;
						export(AxU, oPath + "AxU.dat");
					logger.end();
					if (logger.yes("print result of multiplication")) {
						logger.buf << "A * u = " << AxU;
						logger.log();
					}
					logger.log("check results w/ Mathematica/SymmetricCSlC.nb!");
				logger.end();
				logger.beg("Harwell-Boeing i/o test");
				// …
				logger.end();
			logger.end();
		}
		else if (testNum == 3) {
			/*
				(3) CSR–matrix
			*/
			//logger.beg("CSR matrix test");
			//	ifstream iCSR("Mathematica/matrices/csr.dat"),
			//			 iU("Mathematica/matrices/csr_vector.dat"),
			//			 iV("Mathematica/matrices/csr_vector_t.dat");
			//	ofstream oCSRxU("Mathematica/matrices/csr_mult.dat"),
			//			 oCSRxV("Mathematica/matrices/csr_mult_t.dat");
			//	size_t rows, cols, nnz;
			//	logger.log("loading matrix");
			//	iCSR >> rows >> cols >> nnz;
			//	CSRMatrix<double> CSR(rows, cols, nnz);
			//	iCSR >> CSR;
			//	CSR.save(logger.buf << "fullform:\n");
			//	logger.log();
			//	vector<double> u(cols), v(rows);
			//	iU >> u;
			//	iV >> v;
			//	logger.beg("mult() / multTranspose testing\nresults: Mathematica/matrices/csr_mult(_t).dat");
			//		oCSRxU << CSR * u;
			//		oCSRxV << CSR.t() * v;
			//	logger.end();
			//logger.end();
		}
		else if (testNum == 4) {
			logger.beg("DenseMatrix test");
				logger.beg("create empty matrix");
					Index rows, cols;
					logger.inp("numb of rows", rows)
					      .inp("numb of cols", cols);
					DenseMatrix<double> A(rows, cols);
					logger.buf << "created matrix:\n";
					A.export(logger.buf);
					logger.log();
				logger.end();
				logger.beg("create matrix from ini-list");
					DenseMatrix<double> B {
						{ 1., 2., 7. },
						{ 4., 8., 3. },
						{ 5., 6., 2. },
						{ 3., 2., 1. }
					};
					logger.buf << "created matrix:\n";
					B.export(logger.buf);
					logger.log();
					logger.beg("mult() / multByTranspose() test");
						vector<double> u { 1., 2., 3. }, v { 1., 2., 3., 4. };
						logger.buf << "u       = " << u << '\n'
						           << "v       = " << v << '\n'
						           << "B * u   = " << B * u << '\n'
							       << "B^T * v = " << B.t() * v;
						logger.log();
					logger.end();
					logger.log("set matrix = 3.");
					B = 3.;
					logger.buf << "modified matrix:\n";
					B.export(logger.buf);
					logger.log();
					logger.beg("Gauss elimination w/ partial pivoting");
						B = {
							{ 1.,  4.,  5. },
							{ 7.,  6.,  3. },
							{ 10., 3.5, 1. }
						};
						vector<double> x { 1., 5., 2.3 }, b = B * x;
						logger.buf << "B = \n";
						B.export(logger.buf);
						logger.buf << "b = " << b << '\n'
						           << "x = " << B.GaussElimination(b);
						logger.log();
					logger.end();
				logger.end();
			logger.end();
		}
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}