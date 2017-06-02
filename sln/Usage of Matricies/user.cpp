﻿#include <complex> // for complex–valued matrices
#include "CSCMatrix.hpp"
#include "SymmetricCSlCMatrix.hpp"
#include "SymmetricMatrix.hpp"
#include "DenseMatrix.hpp"
#include "BlockMatrix.hpp"
#include "SymmetricBlockMatrix.hpp"
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
			"DenseMatrix",
			"BlockMatrix",
			"SymmetricBlockMatrix"
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
				logger.beg("precond test");
					double w = 1.5;
					std::vector<double> x(A.getOrder());
					for (Index i = 0; i < x.size(); ++i) x[i] = i + 1.;
					logger.log("x := < 1, 2, ..., n >, w := 1.5");
					logger.log("A := L + D + U, U = L^T");
					logger.buf << "forw subst, [ L + w D]^-1 . x = " << A.forwSubst(x, w) << '\n'
					           << "back subst, [ U + w D]^-1 . x = " << A.backSubst(x, w) << '\n'
					           << "diag subst, [       D]^-1 . x = " << A.diagSubst(x);
					logger.log();
					export(A.forwSubst(x, w), oPath + "forw.dat");
					export(A.backSubst(x, w), oPath + "back.dat");
					export(A.diagSubst(x),    oPath + "diag.dat");
				logger.end();
				logger.beg("Harwell-Boeing i/o test");
				// …
				logger.end();
			logger.end();
		}
		else if (testNum == 3) {
			logger.beg("SymmetricMatrix test");
				logger.beg("create form ini list");
					SymmetricMatrix<double> A {
						11., 12., 13.,
						     22., 23.,
						          33.
					};
					A.export(logger.buf << "created matrix:\n");
					logger.log();
					logger.log("divide by 4.");
					A /= 4.;
					A.export(logger.buf << "modified matrix:\n");
					logger.log();
				logger.end();
			logger.end();
		}
		else if (testNum == 4) {
			logger.beg("DenseMatrix test");
				{
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
					}
					logger.beg("LU decomposition");
						// example from https://people.cs.kuleuven.be/~karl.meerbergen/didactiek/h03g1a/ilu.pdf, p. 18
						DenseMatrix<double> A {
							{	3,	0,	-1,	-1,	0,	-1	},
							{	0,	2,	0,	-1,	0,	0	},
							{	-1,	0,	3,	0,	-1,	0	},
							{	-1,	-1,	0,	2,	0,	-1	},
							{	0,	0,	-1,	0,	3,	-1	},
							{	-1,	0,	0,	-1,	-1,	4	}
						};
						auto b = A * std::vector<double>(6, 1.);
						logger.beg("IKJ / Crout decomposition");
							A.decompositionStrategy = DenseDecompositionStrategy::Crout_L_plus_I_times_D_times_I_plus_U<double>;
							A.decompose().export(logger.buf);
							logger.buf << "\nsoln: " << A.backSubst(A.forwSubst(b, 0.));
							logger.buf << "\nsoln: " << A.backSubst(A.forwSubst(b));
							logger.buf << "\nsoln: " << A.backSubst(A.diagSubst(A.forwSubst(b, 0.)), 0.);
							logger.log();
						logger.end();
						logger.beg("Crout for CSlC");
							CSlCMatrix<double> B {
								{ 0, 3, 4, 7, 9, 10, 10 },
								{ 2, 3, 5, 3, 3, 4, 5, 4, 5, 5 },
								{ -1, -1, -1, -1, 0, -1, 0, 0, -1, -1},
								{ -1, -1, -1, -1, 0, -1, 0, 0, -1, -1 },
								{ 3, 2, 3, 2, 3, 4 }
							};
							B.export(logger.buf);
							logger.log();
							B.decompose().export(logger.buf);
							logger.buf << "\nsoln: " << B.backSubst(B.forwSubst(b, 0.));
							logger.buf << "\nsoln: " << B.backSubst(B.forwSubst(b));
							logger.buf << "\nsoln: " << B.backSubst(B.diagSubst(B.forwSubst(b, 0.)), 0.);
							logger.log();
						logger.end();
						logger.beg("Crout for SymmetricCSlC");
							SymmetricCSlCMatrix<double> S {
								{ 0, 3, 4, 7, 9, 10, 10 },
								{ 2, 3, 5, 3, 3, 4, 5, 4, 5, 5 },
								{ -1, -1, -1, -1, 0, -1, 0, 0, -1, -1},
								{ 3, 2, 3, 2, 3, 4 }
							};
							S.export(logger.buf);
							logger.log();
							S.decompose().export(logger.buf);
							logger.buf << "\nsoln: " << S.backSubst(S.forwSubst(b, 0.));
							logger.buf << "\nsoln: " << S.backSubst(S.forwSubst(b));
							logger.buf << "\nsoln: " << S.backSubst(S.diagSubst(S.forwSubst(b, 0.)), 0.);
							logger.log();
						logger.end();
					logger.end();
				logger.end();
			logger.end();
		} 
		else if (testNum == 5) {
			logger.beg("BlockMatrix test");
				DenseMatrix<int> a {
					{ 1, 2, 3 },
					{ 4, 5, 6 }
				};
				DenseMatrix<int> b {
					{ 1, 2, 3, 6 },
					{ 4, 5, 6, 6 }
				};
				a.export(logger.buf);
				b.export(logger.buf);
				logger.log();
				BlockMatrix<int> A {
					{ nullptr, nullptr },
					{ nullptr, &b },
					{ &a, &b },
				};
				vector<int> u(A.numbOfCols(), 1);
				logger.buf << "block mult:\n" << A * u;
				logger.log();
			logger.end();
		}
		else if (testNum == 6) {
			logger.beg("SymmetricBlockMatrix test");
				// little example of typical saddle point matrix
				// resulting for CFD problems 
				DenseMatrix<int> A11 {
					{ 1, 2, 3 },
					{ 4, 5, 6 },
					{ 7, 8, 9 }
				};
				auto& A22 = A11;
				CSCMatrix<int> A31 {
					{ 0, 1, 1, 3 }, // colptr
					{ 0, 0, 1 }, // rowind
					{ 1, 2, 3 }, // values
					2   // h
				}, A32 {
					{ 0, 0, 1, 3 }, // colptr
					{ 1, 0, 1 }, // rowind
					{ 1, 2, 3 }, // values
					2   // h
				};
				A11.export(logger.buf << "A11 = A22:\n");
				A31.export(logger.buf << "A31:\n");
				A32.export(logger.buf << "A32:\n");
				logger.log();
				SymmetricBlockMatrix<int> A {
					{ &A11, &A22, nullptr }, // diag blocks
					{ // lval blocks
						nullptr, 
						&A31, &A32
					}
				};
				vector<int> u(A.getOrder(), 1);
				logger.buf << "block mult:\n" << A * u;
				logger.log();
			logger.end();
		}
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}