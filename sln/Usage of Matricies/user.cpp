/*
	Alexander Žilykov, Aug 2016
*/

#include <fstream>
#include <string>
#include <complex>
#include "CSRMatrix.hpp"
#include "CSCMatrix.hpp"

int main() {
	// common vars
	string iPath, oPath;
	vector<string> HarwellBoeingMatricies = {
		"mathematica.rra (real    rectangular assembled) -- simple artificial example",
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
						size_t rows, cols, nnz;
						iA >> rows >> cols >> nnz;
						CSCMatrix<double> A(rows, cols, nnz);
						iA >> A;
						vector<double> u(cols), v(rows);
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
						vector<double> AxU(rows), AxV(cols);
						oAxU << (AxU = A * u);
						oAxV << (AxV = A.t() * v);
					logger.end();
					if (logger.yes("print results of multiplications")) {
						logger.buf << "A.u = " << AxU << '\n' << "A^T.v = " << AxV << '\n';
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
					size_t matrixType = logger.opt("choose T for CSCMatrix<T> B", { "double", "complex<double>" }),
					       iMatrixNum = logger.opt("choose input matrix", HarwellBoeingMatricies);
					// leave only *.rra, *.rua etc.
					for_each(HarwellBoeingMatricies.begin(), HarwellBoeingMatricies.end(), [](string& s) { s.resize(s.find_first_of(' ')); });
					string iHBPath = iPath + HarwellBoeingMatricies[iMatrixNum],
					       oHBPath = oPath + HarwellBoeingMatricies[iMatrixNum];
					logger.beg("load info from " + iHBPath);
						logger.beg("load HB header");
							HarwellBoeingHeader header;
							loadHarwellBoeingHeader_f90(iHBPath.c_str(), &header);
							logger.buf << header;
							logger.log();
						logger.end();
						if (matrixType == 0) { // real matrix
							logger.beg("load HB structure");
								CSCMatrix<double> B(header.nrow, header.ncol, header.nnzero);
								B.loadHarwellBoeing(header, iHBPath);
							logger.end();
					logger.end();
							// output
							logger.beg("save matrix to " + oHBPath);
								B.saveHarwellBoeing(oHBPath, { {"title", "My real HB Matrix"} });
							logger.end();
						}
						else { // complex matrix
							logger.beg("load HB structure");
								CSCMatrix<complex<double>> B(header.nrow, header.ncol, header.nnzero);
								B.loadHarwellBoeing(header, iHBPath);
							logger.end();
					logger.end();
							// output
							logger.beg("save matrix to " + oHBPath);
								B.saveHarwellBoeing(oHBPath, { { "title", "My complex HB Matrix" } });
							logger.end();
						}
			logger.end();
		}
		else if (testNum == 3) {
			/*
				(3) CRS–matrix
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