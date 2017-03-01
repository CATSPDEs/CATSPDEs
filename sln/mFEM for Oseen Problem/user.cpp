#include <random> // for initial guess
#include "SingletonLogger.hpp"
#include "MixedFEM.hpp"
#include "constants.hpp"
#include "SymmetricBlockMatrix.hpp" // for saddle point matrix
#include "ProjectionSolvers.hpp"
// Taylor–Hood
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
// Crouzeix–Raviart
#include "Triangle_P1_CrouzeixRaviart.hpp"
#include "Triangle_P0_Lagrange.hpp"

using namespace FEM::Mixed;
using namespace ProjectionSolvers;
using std::string;
using std::vector;



#include "DenseMatrix.hpp"



auto shuffle(Index n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-100., 100.);
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
				diffusion = .01; // Laplace term (inverse Reynolds number)
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
		logger.end();
		auto FEPairIndex = logger.opt("choose FE pair", { "Taylor-Hood", "Crouzeix-Raviart" });
		// Taylor—Hood family
		TriangularScalarFiniteElement
			&velocityFE = Triangle_P2_Lagrange::instance(), 
			&pressureFE = Triangle_P1_Lagrange::instance();
		if (FEPairIndex == 1) {
			// Crouzeix—Raviart family
			velocityFE = Triangle_P1_CrouzeixRaviart::instance();
			pressureFE = Triangle_P0_Lagrange::instance();
		};
		Index numbOfMeshLevels;
		logger.inp("numb of mesh levels (refinements)", numbOfMeshLevels);
		logger.beg("set iterative solver data");
			double eps;
			logger.inp("set eps for solver", eps);
			auto stop = (StoppingCriterion)logger.opt("stopping criterion", { "absolute", "relative" });
			Index i_log;
			logger.inp("log every nth iteration, n", i_log);
		logger.end();
		vector<double> x; // soln vector
		logger.beg("refine mesh");
			Omega.refine(numbOfMeshLevels);
		logger.end();
		logger.beg("assemble system");
			auto system = assembleSystem(PDE, Omega, NeumannBC, DirichletBC, velocityFE, pressureFE);
			auto& A11 = boost::get<0>(system), &A22 = A11;
			auto& B1  = boost::get<1>(system), &B2  = boost::get<2>(system);
			auto& b   = boost::get<3>(system);
			SymmetricBlockMatrix<double> SaddlePointMatrix {
				{ &A11, &A22, nullptr }, // diag
				{ // lval
					nullptr,
					&B1, &B2
				}
			};
		logger.end();
		logger.beg("solve");

				//CSCMatrix<double> SPM;
				//SPM.importHarwellBoeing(oPath + "system/SPM.rua");
				//DenseMatrix<double> SPMd{ SPM };
				//import(b, oPath + "system/b.dat");

			x = Krylov::BiCGStab(SaddlePointMatrix, b, boost::none, eps, stop, 200, 1);
			//x = Krylov::CG(SaddlePointMatrix, b/*, shuffle(SaddlePointMatrix.getOrder()), eps, stop, i_log*/);
			//x = Krylov::BiCGStab(SaddlePointMatrix, b);
			logger.buf << norm(b - SaddlePointMatrix * x);
			logger.log();
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