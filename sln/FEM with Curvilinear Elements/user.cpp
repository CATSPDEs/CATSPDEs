#include "DivGradFEM.hpp"
#include "ProjectionSolvers.hpp"
#include "Triangle_P1_Lagrange.hpp"
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P3_Lagrange.hpp"
#include "FEInterpolant.hpp"
using std::vector;
using std::array;
using std::string;
using std::to_string;
using boost::get;
using namespace FEM::DivGrad;
using namespace ProjectionSolvers;

int main() {
	string iPath("Mathematica/preprocessing/"), oPath("Mathematica/postprocessing/");
	auto& logger = SingletonLogger::instance();
	try {
		logger.beg("build curvilinear mesh");
			Curve2D parabola = [](double const & t) -> Node2D {
				return { 1. - 2. * t, 4. * (1. - t) * t };
			};
			vector<Node2D> nodes {
				{ -1., 0. }, { 0., 0.}, { 0., 1.}, { 1., 0.}
			};
			vector<array<Index, 3>> elements {
				{ 0, 1, 2 }, { 1, 3, 2 }
			};
			vector<array<SignedIndex, 3>> neighbors {
				{ 1, -2, -1 }, { -3, 0, -1 }
			};
			vector<Curve2D> curves { parabola };
			vector<CurvilinearEdge> edges { 
				{ .5, 1., 0 }, { 0., .5, 0 }
			};
			Triangulation Omega {
				nodes, elements, neighbors, curves, edges
			};
			Omega.enumerateRibs();
			Omega.refine(0);
			auto Omega2 = Omega;
			Omega2.deleteCurvilinearEdges();
		logger.end();
		
		logger.beg("Setup pde");
			auto diffusion = [](Node2D const & p) { return 1; };//-32. * ((-1. + p[0])*p[0] + (-1. + p[1])*p[1]); };
			auto reaction = [](Node2D const & p) { return 0; };
			auto force = [](Node2D const & p) { return 2; };
			auto DirichletCondition = [](Node2D const & p) { return 1 - p[0] * p[0] - p[1]; };
			auto strongBCsPredicate = [](Node2D const & p) { return true; };
			auto RobinCoefficient = [](Node2D const & p) { return 0; };
			auto NeumannValue = [](Node2D const & p) { return 0; };
			auto naturalBCsPredicate = [](Node2D const & p) { return false; };
			DiffusionReactionEqn2D PDE { diffusion, reaction, force };
			ScalarBoundaryCondition2D RobinBC{
				{ RobinCoefficient, NeumannValue }, // n . (a ∇u) + R u = N
				naturalBCsPredicate
			}, DirichletBC{ DirichletCondition, strongBCsPredicate };
		logger.end();

		logger.beg("choose FEs");
			vector<TriangularScalarFiniteElement*> FEs {
				&Triangle_P1_Lagrange::instance(),
				&Triangle_P2_Lagrange::instance(),
				&Triangle_P3_Lagrange::instance(),
			};
			auto FEIndex = logger.opt("choose finite element", { "Lagrange P1", "Lagrange P2", "Lagrange P3" });
			auto T_FEIndex = logger.opt("choose finite element for T", { "Lagrange P1", "Lagrange P2", "Lagrange P3" });
			auto& FE = *FEs[FEIndex];
			auto& T_FE = *FEs[T_FEIndex];
		logger.end();

		logger.beg("Assemble matrix");
			auto system = assembleSystem(PDE, Omega, RobinBC, DirichletBC, FE, T_FE);
			auto& A = get<0>(system);
			auto& b = get<1>(system);
		logger.end();
		logger.beg("solve");
			auto maxNumbOfIterations = 1000;
			auto eps = 10e-15;
			auto stop = StoppingCriterion::relative;
			auto i_log = 10;
			auto x = Krylov::CG(A, b, boost::none, maxNumbOfIterations, eps, stop, i_log, 15);
		logger.end();

			auto xInterp = TriangularScalarFEInterpolant(x, FE, Omega);

			logger.buf << xInterp({ -.5,.75 }, 0) << ' ' << xInterp({ -.5,.5 }, 0);
			logger.log();

			Index t;
			auto trueVals = [&](Node2D const & p) {
				return xInterp(p, t);
			};
			auto xInterp2 = TriangularScalarFEInterpolant(trueVals, FE, Omega2, t);


		logger.beg("export soln vector");
			export(x, oPath + "x.dat");
			export(xInterp2.DOFs(), oPath + "y.dat");
		logger.end();
		logger.beg("export");
			Omega.export(oPath + "mesh.ntr", { { "format", "NTR" } });
		logger.end();
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}