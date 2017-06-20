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
// MINI
#include "Triangle_Pt3_LagrangeBubble.hpp"
#include "Triangle_P1_Lagrange.hpp"

using namespace FEM;
using namespace ProjectionSolvers;
using std::string;
using std::vector;
using boost::get;

#include <random> // for initial guess
auto shuffle(Index n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1., 1.);
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
			auto solnIndex = logger.opt("analytic soln", {
				"linear velocity, const pressure",
				// p = 0
				// u = { y, -x }
				"quadratic velocity, const pressure",
				// p = 0
				// u = { -y^2, 2 x y }
				"linear velocity, sin pressure",
				// p = Sin[Pi (x + y)]
				// u = { y, -x }
				"quadratic velocity, sin pressure",
				// p = Sin[Pi (x + y)]
				// u = { -y^2, 2 x y }
				"cubic velocity, quadratic pressure",
				// p = 1/12. + x (1 - x - y)
				// u = { -4 (2 y - 1) (1 - x) x, 4 (2 x - 1) (1 - y) y }
				"non-polynomial velocity and pressure",
				// p = Cos[p[0] p[1]] - mu
				// u = {-p[0] Sin[p[0] p[1]], p[1] Sin[p[0] p[1]]}
				"lid-driven cavity"
				// no analytics
			});
			// BCs			
			auto BCsIndex = logger.opt("choose BCs type", { "only no-slip", "only natural", "mixed" });
			Predicate2D naturalBCPredicate;
			switch (BCsIndex) {
				case 0:
					naturalBCPredicate = [](Node2D const &) { return false; };
					break;
				case 1:
					naturalBCPredicate = [](Node2D const &) { return true; };
					break;
				case 2:
					naturalBCPredicate = [](Node2D const & p) { return p[0] == 1.; };
					break;
			}
			auto problemIndex = logger.opt("choose problem", {
				"Stokes problem",
				"Oseen problem"
			});
			double reaction, diffusion;
			VectorField2D convection, force, NeumannValue, DirichletCondition;
			if (problemIndex == 0) {
				// Stokes
				reaction = 0.; // zero mass term
				diffusion = 1.; // Laplace term
				convection = [](Node2D const &) -> Node2D { return { 0., 0. }; }; // zero wind field
				switch (solnIndex) {
					case 0: // linear velocity, const pressure
						force = [](Node2D const & p) -> Node2D { return { 0., 0. }; };
						NeumannValue = [](Node2D const & p) -> Node2D {
							if (p[0] == 1.) return { 0., -1. };
							if (p[0] == 0.) return { 0., 1. };
							if (p[1] == 1.) return { 1., 0. };
							if (p[1] == 0.) return { -1., 0. };
						};
						DirichletCondition = [](Node2D const & p) -> Node2D { return { p[1], -p[0] }; };
						break;
					case 1: // quadratic velocity, const pressure
						force = [](Node2D const & p) -> Node2D { return { 2., 0. }; };
						NeumannValue = [](Node2D const & p) -> Node2D {
							if (p[0] == 1.) return { -2., 2. * p[1] };
						};
						DirichletCondition = [](Node2D const & p) -> Node2D { return { -p[0] * p[0], 2 * p[0] * p[1] }; };
						break;
					case 2: // linear velocity, sin pressure
						force = [](Node2D const & p) -> Node2D { return { PI * cos(PI * (p[0] + p[1])), PI * cos(PI * (p[0] + p[1])) }; };
						// TODO: add Neumann
						DirichletCondition = [](Node2D const & p) -> Node2D { return { p[1], -p[0] }; };
						break;
					case 3: // quadratic velocity, sin pressure
						force = [](Node2D const & p) -> Node2D { return { 2 + PI * cos(PI * (p[0] + p[1])), PI * cos(PI * (p[0] + p[1])) }; };
						// TODO: add Neumann
						DirichletCondition = [](Node2D const & p) -> Node2D { return { -p[0] * p[0], 2 * p[0] * p[1] }; };
						break;
					case 4: // cubic velocity, quadratic pressure
						force = [](Node2D const & p) -> Node2D { return { 9 - 2 * p[0] - 17 * p[1], -8 + 15 * p[0] }; };
						NeumannValue = [](Node2D const & p) -> Node2D {
							if (p[0] == 1.) return { -49 / 12. + 9 * p[1], -8 * (-1 + p[1]) * p[1] };
							if (p[0] == 0.) return { -47 / 12. + 8 * p[1],  8 * (-1 + p[1]) * p[1] };
							if (p[1] == 1.) return { 8 * (-1 + p[0]) * p[0], 47 / 12. + (-8 + p[0]) * p[0] };
							if (p[1] == 0.) return { -8 * (-1 + p[0]) * p[0], 49 / 12. - p[0] * (7 + p[0]) };
						};
						DirichletCondition = [](Node2D const & p) -> Node2D { return { -4. * (1. - p[0]) * p[0] * (-1. + 2. * p[1]), 4. * (-1. + 2. * p[0]) * (1. - p[1]) * p[1] }; };
						break;
					case 5: // non-polynomial velocity and pressure
						force = [](Node2D const & p) -> Node2D { 
							return { 
								2 * cos(p[0] * p[1]) * p[1] - (pow(p[0], 3) + p[1] + p[0] * p[1] * p[1]) * sin(p[0] * p[1]), 
								-2 * cos(p[0] * p[1]) * p[0] + (-p[0] + p[0] * p[0] * p[1] + pow(p[1], 3)) * sin(p[0] * p[1]) 
							}; 
						};
						NeumannValue = [](Node2D const & p) -> Node2D {
							auto mu = .94608307036718301;
							if (p[0] == 1.) return	{ -cos(p[1]) * (1 + p[1]) - sin(p[1]) + mu, cos(p[1]) * p[1] * p[1] };
							if (p[0] == 0.) return	{ 1 - mu, -p[1] * p[1] };
							if (p[1] == 1.) return	{ -cos(p[0]) * p[0] * p[0], cos(p[0]) * (-1 + p[0]) + sin(p[0]) + mu };
							if (p[1] == 0.) return	{ p[0] * p[0], 1 - mu };
						};
						DirichletCondition = [](Node2D const & p) -> Node2D { return { -p[0] * sin(p[0] * p[1]), p[1] * sin(p[0] * p[1]) }; };
						break;
					case 6: // lid-driven cavity
						// benchmark example from Volker’s “Finite Element Methods for Incompressible Flow Problems,” p. 753
						naturalBCPredicate = [](Node2D const &) { return false; };
						if (BCsIndex) {
							logger.wrn("no-slip BCs will be enforced ");
							BCsIndex = 0;
						}
						force = NeumannValue = [](Node2D const & p) -> Node2D { return { 0., 0. }; };
						auto u1 = [](double x) {
							auto x1 = .1;
							if (x <= x1) return 1. - .25 * (1. - cos(PI * (x1 - x) / x1));
							if (x >= 1. - x1) return 1. - .25 * (1. - cos(PI * (x - (1 - x1)) / x1));
							return 1.;
						};
						DirichletCondition = [=](Node2D const & p) -> Node2D {
							if (p[1] == 1.) return { u1(p[0]), 0. };
							return { 0., 0. }; 
						};
						break;
				}
			}			
			else if (problemIndex == 1) {
				// Oseen problem (same analytics as in Stokes)
				//reaction = 1.; // mass term
				//diffusion = .1; // Laplace term (inverse Reynolds number)
				//convection = [](Node2D const & p) -> Node2D { return { p[0], -p[1] }; }; // wind field
				//force = [&](Node2D const & p) -> Node2D {
				//	return {
				//		 4 * (-2 * p[0] * (-1 + p[1]) + pow(p[0],2)*(-3 + 4 * p[1]) + p[1] * (-1 + p[1] - 4 * diffusion) + 2 * diffusion),
				//		-4 * (-1 + p[1] * (2 + p[1]) + p[0] * (1 - 4 * p[1] - 4 * diffusion) + 2 * diffusion)
				//	};
				//};
				//NeumannValue = [&](Node2D const & p) -> Node2D {
				//	return { 
				//		-4 * (-1 + p[0]) * (-1 + p[1]) * p[1] + 4 * (-1 + 2 * p[0]) * (-1 + 2 * p[1]) * diffusion, 
				//		-8 * (-1 + p[1]) * p[1] * diffusion 
				//	};
				//};
				//DirichletCondition = [](Node2D const & p) -> Node2D {
				//	return {
				//		-4. * (1. - p[0]) * p[0] * (-1. + 2. * p[1]),
				//		 4. * (-1. + 2. * p[0]) * (1. - p[1]) * p[1]
				//	};
				//};
			}
			// PDE
			OseenProblem2D PDE { 
				reaction, convection, diffusion, force, 
				[](Node2D const &) { return 0.; } // div free 
			};
			VectorBoundaryCondition2D 
				NeumannBC	{ NeumannValue, naturalBCPredicate }, 
				DirichletBC	{ DirichletCondition };
			auto FEPairIndex = logger.opt("choose FE pair", { "Taylor-Hood", "Crouzeix-Raviart", "MINI" });
			// Taylor—Hood family
			TriangularScalarFiniteElement
				*velocityFE = &Triangle_P2_Lagrange::instance(),
				*pressureFE = &Triangle_P1_Lagrange::instance();
			// Crouzeix—Raviart family
			if (FEPairIndex == 1) {
				velocityFE = &Triangle_P1_CrouzeixRaviart::instance();
				pressureFE = &Triangle_P0_Lagrange::instance();
			};
			// MINI element
			if (FEPairIndex == 2) {
				velocityFE = &Triangle_Pt3_LagrangeBubble::instance();
				pressureFE = &Triangle_P1_Lagrange::instance();
			};
			vector<string> meshes { "rct_uniform.ntn", "crs_uniform.ntn", "arb_uniform.ntn", "non_uniform.ntn" };
			auto meshesIndex = logger.opt("choose mesh", meshes);
			logger.beg("import initial mesh");
				Triangulation Omega;
				Omega.import(iPath + meshes[meshesIndex]);
				Omega.enumerateRibs();
			logger.end();
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
			auto blockSymmetryType = (BlockSymmetryType)logger.opt("block symmetry type", { "symmetric (indefinite matrix)", "antisymmetric (positive definite matrix)" });		
		logger.end();
		logger.beg("set inner solver data");
			auto precondIndex = logger.opt("choose precond for BiCGStab", { "I", "P_BD", "P_BT" });
			auto Jacobi = Smoothers::relaxedJacobi;
			vector<decltype(Jacobi)> smoothers {
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
			Smoother<CSlCMatrix<double>> smoother = [&](CSlCMatrix<double>& A, vector<double> const & b, vector<double> const & x_0, Index) {
				return smoothers[smoothersIndex](A, b, x_0, omega, nu, 0., StoppingCriterion::absolute, 0);
			};
			auto gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
			++gamma;
			Index numbOfInnerIterations;
			logger.inp("numb of iterations for inner solver", numbOfInnerIterations);
			auto transfer = (TransferType)logger.opt("grid transfer type", { "canonical", "L2" });
		logger.end();
		CSlCMatrix<double> A11;
		CSCMatrix<double> B1, B2;
		vector<double> b;
		logger.beg("build preconditioner");
			logger.beg("(1) MG for Laplace block");
				Multigrid<CSlCMatrix<double>> MG {
					*velocityFE, Omega, numbOfMeshLevels,
					[&](Triangulation const & Omega) {
						auto system = Mixed::assembleSystem(PDE, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE);
						A11 = get<0>(system);
						B1  = get<1>(system);
						B2  = get<2>(system);
						b   = get<3>(system);
						return get<0>(system); // return Laplace block
					},
					transfer
				};
			logger.end();
			logger.beg("(2) pressure mass matrix for Schur complement");
				DiffusionReactionEqn2D ReactionEqn {
					[](Node2D const &) { return 0.; },
					[](Node2D const &) { return 1.; },
					[](Node2D const &) { return 0.; }
				};
				ScalarBoundaryCondition2D 
					rBC { {
							[](Node2D const &) { return 0.; },
							[](Node2D const &) { return 0.; }
						},
						[](Node2D const & p) { return true; }
					}, 
					dBC { [](Node2D const &) { return 0.; } };
				auto pressureMassMatrix = get<0>(
					DivGrad::assembleSystem(ReactionEqn, Omega, rBC, dBC, *pressureFE)
				);
			logger.end();
		logger.end();
		logger.beg("define saddle point matrix");
			Index n = A11.getOrder(), m = B1.numbOfRows();
			if (blockSymmetryType == BlockSymmetryType::antisymmetric) for (Index i = 2 * n; i < b.size(); ++i) b[i] = -b[i];
			SymmetricBlockMatrix<double> SaddlePointMatrix {
				{ &A11, &A11, nullptr }, // diag
				{ // lval
					nullptr,
					&B1, &B2
				},
				blockSymmetryType
			};
		logger.end();
		logger.beg("solve");
			vector<Preconditioner> P {
				/* I  */ [ ](vector<double> const & x) { return x; },
				/* BD */ [&](vector<double> const & x) {
					vector<double>
						x1(x.begin(), x.begin() + n),
						x2(x.begin() + n, x.begin() + 2 * n),
						x3(x.begin() + 2 * n, x.end()),
						y1(n),
						y2(n);
					// (1) Laplace block
					logger.mute = true;
					for (Index i = 0; i < numbOfInnerIterations; ++i) {
						y1 = MG(numbOfMeshLevels, x1, y1, smoother, gamma);
						y2 = MG(numbOfMeshLevels, x2, y2, smoother, gamma);
					}
					logger.mute = false;
					// (2) Schur complement
					auto y3 = pressureMassMatrix.diagSubst(x3);
					if (blockSymmetryType == BlockSymmetryType::antisymmetric) y3 *= -1.;
					// final vector
					vector<double> y;
					y.reserve(2 * n + m);
					y.insert(y.end(), y1.begin(), y1.end());
					y.insert(y.end(), y2.begin(), y2.end());
					y.insert(y.end(), y3.begin(), y3.end());
					return y;
				},
				/* BT */ [&](vector<double> const & x) {
					vector<double>
						x1(x.begin(), x.begin() + n),
						x2(x.begin() + n, x.begin() + 2 * n),
						x3(x.begin() + 2 * n, x.end()),
						y1(n), y2(n);
					// (1) Schur complement
					auto y3 = pressureMassMatrix.diagSubst(x3);
					if (blockSymmetryType == BlockSymmetryType::antisymmetric) y3 *= -1.;
					// (1) Laplace block
					auto f1 = x1 - B1.t() * y3;
					auto f2 = x2 - B2.t() * y3;
					logger.mute = true;
					for (Index i = 0; i < numbOfInnerIterations; ++i) {
						y1 = MG(numbOfMeshLevels, f1, y1, smoother, gamma);
						y2 = MG(numbOfMeshLevels, f2, y2, smoother, gamma);
					}
					logger.mute = false;
					// final vector
					vector<double> y;
					y.reserve(2 * n + m);
					y.insert(y.end(), y1.begin(), y1.end());
					y.insert(y.end(), y2.begin(), y2.end());
					y.insert(y.end(), y3.begin(), y3.end());
					return y;
				}
			};
			auto x = Krylov::PBiCGStab(
				P[precondIndex],
				SaddlePointMatrix, b, 
				boost::none,
				maxNumbOfIterations, eps, stop, i_log
			);
		logger.end();
		logger.beg("export soln vector");
			export(x, oPath + "x.dat");
		logger.end(); 
		logger.beg("export mesh");
			Omega.export(oPath + "mesh.ntr", { { "format", "NTR" } });
		logger.end();
		if (logger.yes("export matrix blocks")) {
			logger.beg("export blocks of matrix and rhs vector");
				static_cast<CSCMatrix<double>>(A11).exportHarwellBoeing(oPath + "system/A11.rua");
				B1.exportHarwellBoeing(oPath + "system/B1.rra");
				B2.exportHarwellBoeing(oPath + "system/B2.rra");
				export(b, oPath + "system/b.dat");
			logger.end();
		}
		logger.exp("stdin.txt");
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}