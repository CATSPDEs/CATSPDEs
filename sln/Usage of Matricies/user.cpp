#include <fstream>
#include <complex>
#include "CSCMatrix.hpp"
#include "SymmetricCSlCMatrix.hpp"
#include "CSRMatrix.hpp"
#include "SingletonLogger.hpp"

/*
	Alexander Žilykov, Aug 2016
*/

int main() {
	// common vars
	string iPath, oPath;
	vector<string> HarwellBoeingMatricies = {
		"mathematica.rra (real    rectangular assembled) -- simple artificial example",
		"mathematica.rsa (real    symmetric   assembled) -- simple artificial example",
		"illc1033.rra    (real    rectangular assembled)",
		"e40r5000.rua    (real    unsymmetric assembled)",
		"qc2534.cua      (comlex  unsymmetric assembled)",
		"young1c.csa     (complex symmetric   assembled)",
		"cegb2802.pse    (pattern symmetric   elemental)"
	};
	// logger
	SingletonLogger& logger = SingletonLogger::instance();
	try {
		size_t testNum = logger.opt(
			"choose matrix format from CATSPDEs collection", 
			{ "CSC", "CSlC", "SymmetricCSlC", "CSR" }
		);
		if (testNum == 0) {
			/*
				(0) CSC–matrix
			*/
			logger.beg("CSC matrix test");
				/*
					(0.0)
				*/
				logger.beg("mult() / multTranspose() test");
					// i/o file pathes
					iPath = "Mathematica/output/CSC/";
					oPath = "output/CSC/";
					string iPathA = iPath + "A.dat", // rect CSC matrix in CATSPDEs sparse format
					       iPathU = iPath + "u.dat", // vector for mult() test
					       iPathV = iPath + "v.dat", // vector for multTranspose() test
					       oPathAxU = oPath + "AxU.dat", 
					       oPathAxV = oPath + "AxV.dat";
					// i/o file streams
					ifstream iA(iPathA),
							 iU(iPathU), 
							 iV(iPathV); 
					ofstream oAxU(oPathAxU), // results of multiplication
							 oAxV(oPathAxV);
					logger.beg("read sparse matrix A from " + iPathA + "\nread vectors u and v from " + iPath);
						CSCMatrix<double> A;
						iA >> A;
						vector<double> u(A.numbOfCols()), v(A.numbOfRows());
						iU >> u;
						iV >> v;
					logger.end();
					if (logger.yes("print A in dense, u, and v")) {
						logger.buf << "A in dense format:\n";
						A.save(logger.buf);
						logger.buf << "\nu = " << u << "\nv = " << v;
						logger.log();
					}
					logger.beg("compute multiplications A.u, A^T.v\nsave results in " + oPath);
						vector<double> AxU(A.numbOfRows()), AxV(A.numbOfCols());
						oAxU << (AxU = A * u);
						oAxV << (AxV = A.t() * v);
					logger.end();
					if (logger.yes("print results of multiplications")) {
						logger.buf << "A.u = " << AxU << "\nA^T.v = " << AxV;
						logger.log();
					}
					logger.log("check results w/ Mathematica/CSC.nb!");
				logger.end();
				/*
					(0.1)
				*/
				logger.beg("Harwell-Boeing i/o test");
					iPath = "HarwellBoeing/";
					oPath = "output/CSC/HarwellBoeing/";
					size_t iMatrixNum = logger.opt("choose input matrix", HarwellBoeingMatricies);
					// leave only *.rra, *.rua etc.
					for_each(HarwellBoeingMatricies.begin(), HarwellBoeingMatricies.end(), [](string& s) { s.resize(s.find_first_of(' ')); });
					string iHBPath = iPath + HarwellBoeingMatricies[iMatrixNum],
					       oHBPath = oPath + HarwellBoeingMatricies[iMatrixNum];
					HarwellBoeingHeader header;
					if (logger.opt("choose T for CSCMatrix<T> B", { "double", "complex<double>" }) == 0) { // real
						logger.beg("load info from " + iHBPath);		
							CSCMatrix<double> B;
							B.loadHarwellBoeing(iHBPath, &header);
							logger.buf << "HB header:\n" << header;
							logger.log();
						logger.end();
						logger.beg("save matrix to " + oHBPath);
							B.saveHarwellBoeing(oHBPath, { { "title", "My real HB matrix" } });
						logger.end();
					}
					else { // complex
						logger.beg("load info from " + iHBPath);
							CSCMatrix<complex<double>> B;
							B.loadHarwellBoeing(iHBPath, &header);
							logger.buf << "HB header:\n" << header;
							logger.log();
						logger.end();
						logger.beg("save matrix to " + oHBPath);
							B.saveHarwellBoeing(oHBPath, { { "title", "My real HB matrix" } });
						logger.end();
					}
				logger.end();
			logger.end();
		}
		else if (testNum == 2) {
			/*
				(2) Symmetric CSlC matrix
			*/
			logger.beg("Symmetric CSlC matrix test");
				/*
					(2.0)
				*/
				logger.beg("mult() test");
					// i/o file pathes
					iPath = "Mathematica/output/SymmetricCSlC/";
					oPath = "output/SymmetricCSlC/";
					string iPathA = iPath + "A.dat", // rect CSC matrix in CATSPDEs sparse format
						   iPathU = iPath + "u.dat", // vector for mult() test
						   oPathAxU = oPath + "AxU.dat";
					// i/o file streams
					ifstream iA(iPathA), iU(iPathU);
					ofstream oAxU(oPathAxU); // results of multiplication
					logger.beg("read sparse matrix A from " + iPathA + "\nread vector u from " + iPathU);
						SymmetricCSlCMatrix<double> A;
						iA >> A;
						vector<double> u(A.getOrder());
						iU >> u;
					logger.end();
					if (logger.yes("print A (in dense format) and u")) {
						logger.buf << "A in dense format:\n";
						A.save(logger.buf);
						logger.buf << "\nu = " << u;
						logger.log();
					}
					logger.beg("compute A.u and save results in " + oPathAxU);
						vector<double> AxU(A.getOrder());
						oAxU << (AxU = A * u);
					logger.end();
					if (logger.yes("print result of multiplication")) {
						logger.buf << "A.u = " << AxU;
						logger.log();
					}
					logger.log("check results w/ Mathematica/SymmetricCSlC.nb!");
				logger.end();
				/*
					(2.1)
				*/
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
	}
	catch (exception const & e) {
		logger.err(e.what());
	}
	return 0;
}