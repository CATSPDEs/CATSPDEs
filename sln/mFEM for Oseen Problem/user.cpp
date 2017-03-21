#include <random> // for initial guess
#include "SingletonLogger.hpp"
#include "MixedFEM.hpp"
#include "constants.hpp"
#include "SymmetricBlockMatrix.hpp" // for saddle point matrix
#include "ProjectionSolvers.hpp" // for outer solver
// for physics-based preconditioner
#include "Multigrid.hpp"
#include "DivGradFEM.hpp"
// Taylor–Hood
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
// Crouzeix–Raviart
#include "Triangle_P1_CrouzeixRaviart.hpp"
#include "Triangle_P0_Lagrange.hpp"

using namespace FEM;
using namespace ProjectionSolvers;
using std::string;
using std::vector;
using boost::get;

auto shuffle(Index n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-10., 10.);
	vector<double> x(n);
	for (auto& el : x) el = dis(gen);
	return x;
}

int main() {
	string iPath("Mathematica/preprocessing/"), oPath("Mathematica/postprocessing/");
	// logger
	auto& logger = SingletonLogger::instance();
	try {
		logger.beg("set up problem");
			double reaction, diffusion;
			VectorField2D convection, force, NeumannValue, DirichletCondition;
			auto problemIndex = logger.opt("choose problem", {
				"Stokes problem",
				"Oseen problem"
			});
			if (problemIndex == 0) {
				// Stokes, analytic soln:
				// p = 4 y (1 - x) (1 - y)
				// u = { -4 (2 y - 1) (1 - x) x, 4 (2 x - 1) (1 - y) y }
				reaction = 0.; // zero mass term
				diffusion = 1.; // Laplace term
				convection = [](Node2D const &) -> Node2D { return { 0., 0. }; }; // zero wind field
				force = [](Node2D const & p) -> Node2D { 
					return { 
						 4. * (2. + (p[1] - 5.) * p[1]), 
						-4. - 8. * p[1] + 4. * p[0] * (3. + 2. * p[1])
					};
				};
				NeumannValue = [](Node2D const & p) -> Node2D {
					return {
						-4. * (-1. - (-3. + p[1]) * p[1] + p[0] * (2. + (-5. + p[1]) * p[1])),
						-8. * (-1 + p[1]) * p[1]
					};
				};
				DirichletCondition = [](Node2D const & p) -> Node2D {
					return {
						-4. * (1. - p[0]) * p[0] * (-1. + 2. * p[1]), 
						 4. * (-1. + 2. * p[0]) * (1. - p[1]) * p[1]
					};
				};
			}			
			else if (problemIndex == 1) {
				// Oseen problem (same analytics as in Stokes)
				reaction = 1.; // mass term
				diffusion = .1; // Laplace term (inverse Reynolds number)
				convection = [](Node2D const & p) -> Node2D { return { p[0], -p[1] }; }; // wind field
				force = [&](Node2D const & p) -> Node2D {
					return {
						 4 * (-2 * p[0] * (-1 + p[1]) + pow(p[0],2)*(-3 + 4 * p[1]) + p[1] * (-1 + p[1] - 4 * diffusion) + 2 * diffusion),
						-4 * (-1 + p[1] * (2 + p[1]) + p[0] * (1 - 4 * p[1] - 4 * diffusion) + 2 * diffusion)
					};
				};
				NeumannValue = [&](Node2D const & p) -> Node2D {
					return { 
						-4 * (-1 + p[0]) * (-1 + p[1]) * p[1] + 4 * (-1 + 2 * p[0]) * (-1 + 2 * p[1]) * diffusion, 
						-8 * (-1 + p[1]) * p[1] * diffusion 
					};
				};
				DirichletCondition = [](Node2D const & p) -> Node2D {
					return {
						-4. * (1. - p[0]) * p[0] * (-1. + 2. * p[1]),
						 4. * (-1. + 2. * p[0]) * (1. - p[1]) * p[1]
					};
				};
			}
			// PDE
			OseenProblem2D PDE { 
				reaction, convection, diffusion, force, 
				[](Node2D const &) { return 0.; } // div free 
			};
			// BCs
			VectorBoundaryCondition2D NeumannBC {
				NeumannValue,
				[](Node2D const & p) { return p[0] == 1.; }
			}, DirichletBC { DirichletCondition };
			logger.beg("import initial mesh");
				Triangulation Omega;
				Omega.import(iPath + "mesh.ntn");
				Omega.enumerateRibs();
			logger.end();
			auto FEPairIndex = logger.opt("choose FE pair", { "Taylor-Hood", "Crouzeix-Raviart" });
			// Taylor—Hood family
			TriangularScalarFiniteElement
				*velocityFE = &Triangle_P2_Lagrange::instance(),
				*pressureFE = &Triangle_P1_Lagrange::instance();
			if (FEPairIndex == 1) {
				// Crouzeix—Raviart family
				velocityFE = &Triangle_P1_CrouzeixRaviart::instance();
				pressureFE = &Triangle_P0_Lagrange::instance();
			};
			Index numbOfMeshLevels;
			logger.inp("numb of mesh levels (refinements)", numbOfMeshLevels);
		logger.end();
		logger.beg("set outer solver data");
			double eps;
			logger.inp("set eps for solver", eps);
			auto stop = (StoppingCriterion)logger.opt("stopping criterion", { "absolute", "relative" });
			Index maxNumbOfIterations;
			logger.inp("max numb of iterations", maxNumbOfIterations);
			Index i_log;
			logger.inp("log every nth iteration, n", i_log);
		logger.end();
		logger.beg("set inner solver data (BD-preconditioner)");
			auto innerSolver = Smoothers::relaxedJacobi;
			vector<decltype(innerSolver)> smoothers {
				Smoothers::relaxedJacobi,
				Smoothers::forwSOR,
				Smoothers::backSOR,
				Smoothers::SSOR
			};
			auto smoothersIndex = logger.opt("smoothing technique", { "relaxed Jacobi", "forward SOR", "backward SOR", "SSOR" });
			Index nu;
			logger.inp("numb of pre- and post-smoothing iterations", nu);
			double omega;
			logger.inp("relaxation parameter", omega);
			Smoother<CSlCMatrix<double>> smoother = [&](CSlCMatrix<double>& A, vector<double> const & b, vector<double> const & x_0) {
				return smoothers[smoothersIndex](A, b, x_0, omega, nu, 0., StoppingCriterion::absolute, 0);
			};
			auto gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
			++gamma;
			Index numbOfInnerIterations;
			logger.inp("numb of iterations for inner solver", numbOfInnerIterations);
			auto transfer = (TransferType)logger.opt("grid transfer type", { "canonical", "L2" });
		logger.end();
		logger.beg("build preconditioner");
			Multigrid<CSlCMatrix<double>> MG {
				*velocityFE, Omega, numbOfMeshLevels,
				[&](Triangulation const & Omega) {
					auto system = Mixed::assembleSystem(PDE, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE);
					return get<0>(system); // return Laplace block
				},
				transfer
			};
			DiffusionReactionEqn2D ReactionEqn {
				[](Node2D const &) { return 0.; },
				[](Node2D const &) { return 1.; },
				[](Node2D const &) { return 0.; }
			};
			ScalarBoundaryCondition2D rBC { {
					[](Node2D const &) { return 0.; },
					[](Node2D const &) { return 0.; }
				},
				[](Node2D const & p) { return true; }
			}, dBC {
				[](Node2D const &) { return 0.; }
			};
			auto pressureMassMatrix = get<0>(
				DivGrad::assembleSystem(ReactionEqn, Omega, rBC, dBC, *pressureFE)
			);
		logger.end();
		logger.beg("assemble saddle point system");
			auto system = Mixed::assembleSystem(PDE, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE);
			auto& A11 = get<0>(system), &A22 = A11;
			auto& B1  = get<1>(system), &B2  = get<2>(system);
			auto& b   = get<3>(system);
			SymmetricBlockMatrix<double> SaddlePointMatrix {
				{ &A11, &A22, nullptr }, // diag
				{ // lval
					nullptr,
					&B1, &B2
				}
			};
		logger.end();
		logger.beg("solve");
			Preconditioner BlockDiagonalPreconditioner = [&](vector<double> const & x) {
				vector<double>
					x1(x.begin(), x.begin() + A11.getOrder()),
					x2(x.begin() + A11.getOrder(), x.begin() + 2 * A11.getOrder()),
					x3(x.begin() + 2 * A11.getOrder(), x.end()),
					y1(A11.getOrder()), 
					y2(A11.getOrder()), 
					y3(pressureMassMatrix.getOrder());
				// (1) Laplace block
				logger.mute = true;
				for (Index i = 0; i < numbOfInnerIterations; ++i) {
					y1 = MG(numbOfMeshLevels, x1, y1, smoother, gamma);
					y2 = MG(numbOfMeshLevels, x2, y2, smoother, gamma);
				}
				logger.mute = false;
				// (2) Schur complement
				y3 = pressureMassMatrix.diagSubst(x3);
				// final vector
				vector<double> y;
				y.reserve(SaddlePointMatrix.getOrder());
				y.insert(y.end(), y1.begin(), y1.end());
				y.insert(y.end(), y2.begin(), y2.end());
				y.insert(y.end(), y3.begin(), y3.end());
				return y;
			};
			auto x = Krylov::PBiCGStab(
				BlockDiagonalPreconditioner,
				SaddlePointMatrix, b, shuffle(SaddlePointMatrix.getOrder()), maxNumbOfIterations, eps, stop, i_log
			);
			//CSCMatrix<double> SPM;
			//SPM.importHarwellBoeing(oPath + "system/SPM.rua");
			//DenseMatrix<double> SPMd{ SPM };
			//import(b, oPath + "system/b.dat");
			//x = Krylov::CG(SaddlePointMatrix, b, boost::none, 10000, eps, stop, 10, 0);
			//x = Krylov::CG(SaddlePointMatrix, b/*, shuffle(SaddlePointMatrix.getOrder()), eps, stop, i_log*/);
			//x = Krylov::BiCGStab(SaddlePointMatrix, b);
			//logger.buf << norm(b - SaddlePointMatrix * x);
			//logger.log();
		logger.end();
		logger.beg("export soln vector");
			export(x, oPath + "x.dat");
		logger.end(); 
		logger.beg("export mesh");
			Omega.export(oPath + "mesh.ntr", { { "format", "NTR" } });
		logger.end();
		logger.beg("export blocks of matrix and rhs vector");
			static_cast<CSCMatrix<double>>(A11).exportHarwellBoeing(oPath + "system/A11.rua");
			B1.exportHarwellBoeing(oPath + "system/B1.rra");
			B2.exportHarwellBoeing(oPath + "system/B2.rra");
			export(b, oPath + "system/b.dat");
		logger.end();
		logger.exp("stdin.txt");
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}