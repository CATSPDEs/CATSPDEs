#include "SingletonLogger.hpp"
#include "MixedFEM.hpp"
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"

using namespace FEM::Mixed;
using std::string;

int main() {
	string iPath("Mathematica/");
	// logger
	auto& logger = SingletonLogger::instance();
	try {
		// set PDE
		OseenProblem2D PDE {
			1.,
			[](Node2D const & p) -> Node2D {
				return { p[0], p[1] };
			},
			2.,
			[](Node2D const & p) -> Node2D {
				return { 1., 1. };
			},
			[](Node2D const & p) {
				return norm(p);
			}
		};
		// import and set up mesh
		Triangulation Omega;
		Omega.import(iPath + "mesh.ntn");
		Omega.enumerateRibs()
		     .export(iPath + "mesh.ntr", { {"format", "NTR"} });
		// set up BCs
		VectorBoundaryCondition2D DirichletBC {
			[](Node2D const &) -> Node2D {
				return { 0., 0. };
			}
		}, NeumannBC {
			[](Node2D const &) -> Node2D {
				return { 1., 2. };
			},
			[](Node2D const & p) {
				return p[0] == 10.;
			}
		};
		// choose FEs
		auto& velocityFE = Triangle_P2_Lagrange::instance();
		auto& pressureFE = Triangle_P1_Lagrange::instance();
		// assmeble system
		assembleSystem(PDE, Omega, DirichletBC, NeumannBC, velocityFE, pressureFE);
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}