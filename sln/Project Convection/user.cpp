#include <memory>
// finite elements to use
#include "Triangle_P0_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_CrouzeixRaviart.hpp"
// interpolant for diffusion
#include "FEInterpolant.hpp"
// assembler
#include "DivGradFEM.hpp"
// solvers
#include "ProjectionSolvers.hpp"
// pi
#include "constants.hpp"

int main() {
	string iPath("Mathematica/preprocessing/"), oPath("Mathematica/postprocessing/");
	auto& logger = SingletonLogger::instance();
	try {
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}