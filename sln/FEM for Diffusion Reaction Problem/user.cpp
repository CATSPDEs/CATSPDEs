﻿#include <random> // for initial guess in MG
#include "Multigrid.hpp"

// finite elements to use
#include "Triangle_P1_Lagrange.hpp"
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_CrouzeixRaviart.hpp"

#include "constants.hpp"

using std::vector;
using std::string;
using std::to_string;
using boost::get;
using namespace FEM::DivGrad;
using namespace ProjectionSolvers;

// shuffle initial guess 
// in order to introduce all possible freq in the error, thus test MG
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
	auto& logger = SingletonLogger::instance();
	try {
		logger.beg("set up problem");
			ScalarField2D diffusion, reaction, force, 
			              NeumannValue, RobinCoefficient, DirichletCondition;
			Predicate2D   naturalBCsPredicate, strongBCsPredicate;
			auto problemIndex = logger.opt("choose problem", {
				"Laplace w/ hom. Dirichlet BCs (zero analytical soln)",
				"Poisson w/ inhom. Robin BCs",
				"Diffusion-Reaction w/ mixed Dirichlet and Robin BCs"
			});
			if (problemIndex == 0) {
				// Laplace,
				// analytic soln = 0
				diffusion = [](Node2D const) { return 1.; };
				reaction = force = NeumannValue = RobinCoefficient = DirichletCondition = [](Node2D const) { return 0.; };
				naturalBCsPredicate = [](Node2D const) { return false; };
				strongBCsPredicate  = [](Node2D const) { return true; };
			}
			else if (problemIndex == 1) {
				// Poisson,
				// analytic soln = 16 x y (x - 1) (y - 1)
				diffusion = [](Node2D const) { return 1.; };
				reaction = DirichletCondition = [](Node2D const) { return 0.; };
				force = [](Node2D const & p) { return -32. * ((-1. + p[0])*p[0] + (-1. + p[1])*p[1]); };
				RobinCoefficient = [](Node2D const) { return 1.; };
				NeumannValue = [](Node2D const & p) { 
					if (p[0] == 0. || p[0] == 1.) return 16. * (-1. + p[1]) * p[1];
					if (p[1] == 0. || p[1] == 1.) return 16. * (-1. + p[0]) * p[0];
				};
				naturalBCsPredicate = [](Node2D const) { return true; };
				strongBCsPredicate  = [](Node2D const) { return false; };
			}
			else if (problemIndex == 2) {
				// Diffusion-Reaction,
				// analytic soln = Sin[2 Pi x y]
				diffusion = reaction = RobinCoefficient = [](Node2D const) { return 1.; };
				DirichletCondition = [](Node2D const & p) { return sin(2. * PI * p[0] * p[1]); };
				force = [](Node2D const & p) { return (1. + 4. * PI * PI * p * p) * sin(2. * PI * p[0] * p[1]); };
				NeumannValue = [](Node2D const & p) {
					if (p[0] == 1.) return  2. * PI * cos(2. * PI * p[1]) * p[1] + sin(2. * PI * p[1]);
					if (p[1] == 0.) return -2. * PI * p[0];
				};
				naturalBCsPredicate = [](Node2D const & p) { return p[0] == 1. || p[1] == 0.; }; // _|-part of square
				strongBCsPredicate  = [](Node2D const & p) { return p[0] == 0. || p[1] == 1.; }; // Г-part
			}
			// PDE
			DiffusionReactionEqn2D PDE { diffusion, reaction, force };
			// BCs
			ScalarBoundaryCondition2D RobinBC { 
				{ RobinCoefficient, NeumannValue }, // n . (a ∇u) + R u = N
				naturalBCsPredicate
			}, DirichletBC { DirichletCondition, strongBCsPredicate };
			// mesh
			Index initialRefCount;
			vector<string> meshType { "uniform.ntn", "nonuniform.ntn" };
			auto meshTypeIndex = logger.opt("mesh type", meshType);
			logger.inp("refine initial mesh n times, n", initialRefCount);
			logger.beg("import initial mesh");
				Triangulation Omega;
				Omega.import(iPath + meshType[meshTypeIndex]).refine(initialRefCount);
				logger.beg("enumerate ribs");
					Omega.enumerateRibs();
				logger.end();
			logger.end();
			// FE
			vector<TriangularScalarFiniteElement*> FEs {
				&Triangle_P1_Lagrange::instance(),
				&Triangle_P2_Lagrange::instance(),
				&Triangle_P1_CrouzeixRaviart::instance()
			};
			auto FEIndex = logger.opt("choose finite element", { "Lagrange P1", "Lagrange P2", "Crouzeix-Raviart P1" });
			auto& FE = *FEs[FEIndex];
			// numb of refinements
			Index numbOfMeshLevels;
			logger.inp("numb of mesh levels (refinements)", numbOfMeshLevels);
		logger.end();
		auto method = logger.opt("choose solving technique", {
			"MG",
			"PCG w/ MG as an inner iteration",
			"CG",
			"Smoother"
		});
		logger.beg("set iterative solver data");
			double eps;
			logger.inp("set eps for solver", eps);
			auto stop = (StoppingCriterion)logger.opt("stopping criterion", { "absolute", "relative" });
			Index i_log;
			logger.inp("log every nth iteration, n", i_log);
		logger.end();
		vector<double> x; // soln vector
		if (method == 0) {
			logger.beg("setup multigrid data");
				Index maxNumbOfIterations;
				logger.inp("max numb of iterations", maxNumbOfIterations);
				vector<string> smoothersNames { "relaxed Jacobi", "forward SOR", "backward SOR", "SSOR" };
				auto dummy = Smoothers::relaxedJacobi;
				vector<decltype(dummy)> smoothers {
					Smoothers::relaxedJacobi,
					Smoothers::forwSOR,
					Smoothers::backSOR,
					Smoothers::SSOR
				};
				auto smoothersIndex = logger.opt("smoothing technique", smoothersNames);
				Index nu;
				logger.inp("numb of pre- and post-smoothing iterations", nu);
				double omega;
				logger.inp("relaxation parameter", omega);
				Smoother<SymmetricCSlCMatrix<double>> smoother = [&](SymmetricCSlCMatrix<double>& A, std::vector<double> const & b, std::vector<double> const & x_0) {
					return smoothers[smoothersIndex](A, b, x_0, omega, nu, 0., StoppingCriterion::absolute, 0);
				};
				auto gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
				++gamma;
				auto transfer = (TransferType)logger.opt("grid transfer type", { "canonical", "L2" });
				Multigrid<SymmetricCSlCMatrix<double>> MG {
					FE, Omega, numbOfMeshLevels,
					[&](Triangulation const & Omega) {
						return assembleSystem(PDE, Omega, RobinBC, DirichletBC, FE);
					},
					transfer
				};
				auto& A = MG.A();
				auto& b = MG.b();
			logger.end();
			logger.beg("solve w/ MG");
				x = shuffle(A.getOrder()); // import(x, oPath + "shuffled guesses/x_0.dat");
				// initial residual
				double r_0 = norm(b - A * x), r = r_0, r_new;
				logInitialResidual(r_0, maxNumbOfIterations);
				// MG as a stand-alone iteration
				Index i;
				for (i = 1; i <= maxNumbOfIterations; ++i) {
					logger.mute = true;
					x = MG(numbOfMeshLevels, b, x, smoother, gamma);
					logger.mute = false;
					r_new = norm(b - A * x);
					if (i_log && i % i_log == 0) logResidualReduction(r, r_new, i, maxNumbOfIterations);
					r = r_new;
					if      (stop == StoppingCriterion::absolute && r       < eps) break;
					else if (stop == StoppingCriterion::relative && r / r_0 < eps) break;
				}
				logFinalResidual(r_0, r_new, i > maxNumbOfIterations ? maxNumbOfIterations : i, maxNumbOfIterations);
				if (i > maxNumbOfIterations) logger.wrn("MG failed");
			logger.end();
			if (transfer == TransferType::canonical && logger.yes("export prolongations and system matrices")) {
				logger.beg("export");
					for (Index i = 0; i <= numbOfMeshLevels; ++i)
						static_cast<CSCMatrix<double>>(static_cast<CSlCMatrix<double>>(MG.A(i))).exportHarwellBoeing(oPath + "mg/A" + to_string(i) + ".rsa");
					for (Index i = 0; i < numbOfMeshLevels; ++i)
						MG.P_Hh(i).exportHarwellBoeing(oPath + "mg/P" + to_string(i) + ".rra");
				logger.end();
			}
			// save stdin commands
			logger.exp("mg.txt");
		}
		else if (method == 1) {
			logger.beg("setup multigrid data");
				Index maxNumbOfIterations;
				logger.inp("max numb of iterations", maxNumbOfIterations);
				vector<string> smoothersNames { "relaxed Jacobi", "forward SOR", "backward SOR", "SSOR" };
				auto dummy = Smoothers::relaxedJacobi;
				vector<decltype(dummy)> smoothers{
					Smoothers::relaxedJacobi,
					Smoothers::forwSOR,
					Smoothers::backSOR,
					Smoothers::SSOR
				};
				auto smoothersIndex = logger.opt("smoothing technique", smoothersNames);
				Index nu;
				logger.inp("numb of pre- and post-smoothing iterations", nu);
				double omega;
				logger.inp("relaxation parameter", omega);
				Smoother<SymmetricCSlCMatrix<double>> smoother = [&](SymmetricCSlCMatrix<double>& A, std::vector<double> const & b, std::vector<double> const & x_0) {
					return smoothers[smoothersIndex](A, b, x_0, omega, nu, 0., StoppingCriterion::absolute, 0);
				};
				auto gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
				++gamma;
				Index numbOfInnerIterations;
				logger.inp("numb of iterations for inner solver", numbOfInnerIterations);
				auto transfer = (TransferType)logger.opt("grid transfer type", { "canonical", "L2" });
				Multigrid<SymmetricCSlCMatrix<double>> MG {
					FE, Omega, numbOfMeshLevels,
					[&](Triangulation const & Omega) {
						return assembleSystem(PDE, Omega, RobinBC, DirichletBC, FE);
					},
					transfer
				};
				auto& A = MG.A();
				auto& b = MG.b();
			logger.end();
			logger.beg("solve w/ PCG");
				x = Krylov::PCG([&](std::vector<double> const & x) {
					std::vector<double> y(A.getOrder());
					logger.mute = true;
					for (Index i = 0; i < numbOfInnerIterations; ++i)
						y = MG(numbOfMeshLevels, x, y, smoother, gamma);
					logger.mute = false;
					return y;
				}, A, b, shuffle(A.getOrder()), eps, stop, i_log);
			logger.end();
			// save stdin commands
			logger.exp("pcg.txt");
		}
		else if (method == 2) {
			logger.beg("refine mesh");
				Omega.refine(numbOfMeshLevels);
			logger.end();
			logger.beg("assemble system");
				auto system = assembleSystem(PDE, Omega, RobinBC, DirichletBC, FE);
				auto& A = get<0>(system);
				auto& b = get<1>(system);
			logger.end();
			logger.beg("solve w/ CG");
				x = Krylov::CG(A, b, shuffle(A.getOrder()), eps, stop, i_log);
			logger.end();
			// save stdin commands
			logger.exp("cg.txt");
		}
		else if (method == 3) {
			logger.beg("setup solver data");
				vector<string> smoothersNames { "relaxed Jacobi", "forward SOR", "backward SOR", "SSOR" };
				using namespace Smoothers;
				auto dummy = relaxedJacobi;
				vector<decltype(dummy)> smoothers { relaxedJacobi, forwSOR, backSOR, SSOR };
				auto smoothersIndex = logger.opt("choose smoother", smoothersNames);
				Index numbOfIterations;
				logger.inp("numb of iterations", numbOfIterations);
				double omega;
				logger.inp("relaxation parameter", omega);
			logger.end();
			logger.beg("refine mesh");
				Omega.refine(numbOfMeshLevels);
			logger.end();
			logger.beg("assemble system");
				auto system = assembleSystem(PDE, Omega, RobinBC, DirichletBC, FE);
				auto& A = get<0>(system);
				auto& b = get<1>(system);
			logger.end();
			if (logger.opt("where to log", {"stdout", "save to file"})) {
				logger.beg("smooth");
					x = shuffle(A.getOrder());
					vector<decltype(x)> residuals;
					residuals.reserve(numbOfIterations + 1);
					residuals.emplace_back(A * x - b);
					logger.mute = true;
					for (Index i = 1; i <= numbOfIterations; ++i) {
						x = smoothers[smoothersIndex](A, b, x, omega, 1, eps, stop, 0);
						residuals.emplace_back(b - A * x);
					}
					logger.mute = false;
				logger.end();
				logger.beg("export residuals history");
					export(residuals, oPath + "residuals history/" + smoothersNames[smoothersIndex] + ".dat");
				logger.end();
			}
			else x = smoothers[smoothersIndex](A, b, shuffle(A.getOrder()), omega, numbOfIterations, eps, stop, i_log);
			logger.exp("smoothers.txt");
		}
		logger.beg("export soln vector");
			export(x, oPath + "x.dat");
		logger.end();
		logger.beg("export mesh");
			Omega.export(oPath + "mesh.ntr", { {"format", "NTR"} });
		logger.end();
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}