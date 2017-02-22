#include "SingletonLogger.hpp"
#include "MixedFEM.hpp"
#include "constants.hpp"

#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"

#include "Triangle_P1_CrouzeixRaviart.hpp"
#include "Triangle_P0_Lagrange.hpp"

#include "DenseMatrix.hpp"

using namespace FEM::Mixed;
using std::string;

int main() {
	string iPath("Mathematica/preprocessing/"), oPath("Mathematica/postprocessing/");
	// logger
	auto& logger = SingletonLogger::instance();
	try {
		// set PDE
		OseenProblem2D PDE {
			1., // mass coef
			[](Node2D const & p) -> Node2D { // wind field
				return { p[0], -p[1] };
			},
			.01, // inv Reynolds
			[](Node2D const & p) -> Node2D { // ext force field
				return {
					0.08 - 0.16*p[1] + p[0] * (8. - 8.*p[1] + p[0] * (-12. + 16.*p[1])) - 3.141592653589793*sin(3.141592653589793*(p[0] + p[1])),
					-0.08 - 4.*p[1] *p[1] + p[0] * (0.16 + 8.*p[1]) - 3.141592653589793*sin(3.141592653589793*(p[0] + p[1]))
				};
			},
			[](Node2D const &) { // div free
				return 0.;
			}
		};
		// import and set up mesh
		Triangulation Omega;
		Omega.import(iPath + "mesh.ntn");
		Omega.enumerateRibs().refine(1).export(oPath + "mesh.ntr", { {"format", "NTR"} });
		// set up BCs
		VectorBoundaryCondition2D NeumannBC {
			[](Node2D const & p) -> Node2D {
				return { -0.96 - 1.*cos(3.141592653589793*(p[0] + p[1])) + p[0] * (-0.08 + 0.16*p[1]) - 0.08*p[1], -0.08*(-1. + p[1])*p[1] };
				//if (p[0] == 1.) return { -2. - 16.*p[0] - 16.*p[1] + 32.*p[0] * p[1] - sin(5.*(p[0] + p[1])), -16.*(-1. + p[1])*p[1] };
				//if (p[0] == 0.) return { 2. + p[0] * (16. - 32.*p[1]) + 16.*p[1] + sin(5.*(p[0] + p[1])), 16.*(-1. + p[1])*p[1] };
				//if (p[1] == 1.) return { 16.*(-1. + p[0])*p[0], 2.*(-9. + p[0] * (8. - 16.*p[1]) + 8.*p[1]) - sin(5.*(p[0] + p[1])) };
				//if (p[1] == 0.) return { -16.*(-1. + p[0])*p[0], 2.*(9. - 8.*p[1] + 8.*p[0] * (-1. + 2.*p[1])) + sin(5.*(p[0] + p[1])) };
			},
			[](Node2D const & p) {
				return p[0] == 1.;
			}
		}, DirichletBC {
			[&](Node2D const & p) -> Node2D {
				return { 
					-4.*(1. - 1.*p[0])*p[0] * (-1. + 2.*p[1]), 
					4.*(-1. + 2.*p[0])*(1. - 1.*p[1])*p[1] 
				};
			}
		};
		// choose FEs
		//auto& velocityFE = Triangle_P1_CrouzeixRaviart::instance();
		//auto& pressureFE = Triangle_P0_Lagrange::instance();
		auto& velocityFE = Triangle_P2_Lagrange::instance();
		auto& pressureFE = Triangle_P1_Lagrange::instance();
		// assmeble system
		auto system = assembleSystem(PDE, Omega, NeumannBC, DirichletBC, velocityFE, pressureFE);
		auto& mtx = boost::get<0>(system);

		static_cast<decltype(mtx.B1)>(mtx.A11).exportHarwellBoeing(oPath + "A11.rua");
		mtx.B1.exportHarwellBoeing(oPath + "B1.rra");
		mtx.B2.exportHarwellBoeing(oPath + "B2.rra");
		export(boost::get<1>(system), oPath + "b.dat");

	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}