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
		"illc1033.rra (real    rectangular assembled)",
		"e40r5000.rua (real    unsymmetric assembled)",
		"qc2534.cua   (comlex  unsymmetric assembled)",
		"young1c.csa  (complex symmetric   assembled)",
		"cegb2802.pse (pattern symmetric   elemental)"
	};
	size_t i;
	// logger
	SingletonLogger& logger = SingletonLogger::instance();
	try {
		size_t testNum = logger.opt(
			"choose matrix format from CATSPDEs collection", 
			{ "CSC", "CSlC", "SymmetricCSlC", "CSR" }
		);
		/*
			(0) CSC–matrix
		*/
		if (testNum == 0) {
			logger.beg("CSC matrix test");
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
				logger.beg("Harwell-Boeing i/o test");
					iPath = "HarwellBoeing/";
					size_t iMatrixNum = logger.opt("choose input matrix", HarwellBoeingMatricies);
					// leave only *.rra, *.rua etc.
					for_each(HarwellBoeingMatricies.begin(), HarwellBoeingMatricies.end(), [](string& s) { s.resize(s.find_first_of(' ')); });
					logger.beg("load HB header");
						string iHBPath = iPath + HarwellBoeingMatricies[iMatrixNum];
						HarwellBoeingHeader header;
						loadHarwellBoeingHeader_f90(iHBPath.c_str(), &header);
						logger.buf << "header of " << iHBPath << ":\n" << header;
						logger.log();
					logger.end();
					if (logger.opt("choose T for CSCMatrix<T> B", { "double", "complex<double>" }) == 0) {
						// real
						CSCMatrix<double> B(header);
						B.loadHarwellBoeing(iHBPath);
					}
					else {
						// complex
						CSCMatrix<complex<double>> B(header);
						B.loadHarwellBoeing(iHBPath);
					}
				logger.end();
			logger.end();
		}
		/*
		(3) CRS–matrix
		*/
		else if (testNum == 3) {
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