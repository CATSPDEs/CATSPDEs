#include <cmath> // pow
#include "SingletonLogger.hpp"
// FEs to build prolongation matrices for
#include "Triangle_P0_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_CrouzeixRaviart.hpp"
#include "Triangle_Pt3_LagrangeBubble.hpp"
// interpolant
#include "FEInterpolant.hpp"
// prolongation matrix data structure
#include "CSCMatrix.hpp"

#include "constants.hpp"

using std::string;
using std::vector;

int main() {
	auto& logger = SingletonLogger::instance();
	string path("Mathematica/");
	try {
		logger.beg("choose FE");
			vector<TriangularScalarFiniteElement*> FEs {
				&Triangle_P0_Lagrange::instance(),
				&Triangle_P1_Lagrange::instance(),
				&Triangle_P2_Lagrange::instance(),
				&Triangle_P1_CrouzeixRaviart::instance(),
				&Triangle_Pt3_LagrangeBubble::instance()
			};
			auto FEIndex = logger.opt("choose finite element", { "Lagrange P0", "Lagrange P1", "Lagrange P2", "Crouzeix-Raviart", "Lagrange P1 Bubble" });
			auto& FE = *FEs[FEIndex];
		logger.end();
		logger.beg("import mesh");
			Triangulation Omega;
			Omega.import(path + "arb_uniform.ntn");
			//Omega.refine();
			Omega.enumerateRibs();
			Omega.export(path + "mesh0.ntr", { { "format", "NTR" } });
		logger.end();
		logger.beg("create interp");
			vector<ScalarField2D> funcs {
				[](Node2D const & p) { return sin(2. * PI * (1. - p[0]) * (1. - p[1])); },
				[](Node2D const & p) { return 50. * p[0] * p[1] * pow(p[0] - 1., 2.) * pow(1. - p[1], 2.); }
			};
			auto& func = funcs[logger.opt("choose func", { "sine", "dome" })];
			TriangularScalarFEInterpolant interp { func, FE, Omega };
			export(interp.DOFs(), path + "x.dat");
		logger.end();
		logger.beg("non-uniform refinement");
			Index n;
			logger.inp("numb of elements to refine", n);
			vector<Index> elems2refine_untreated(n);
			logger.inp("indicies of elements to refine", elems2refine_untreated);
			Indicies elems2refine(elems2refine_untreated.begin(), elems2refine_untreated.end());
			auto OmegaCoarse = Omega; // copy
			Omega.refine(elems2refine);
			auto P0 = FE.prolongate(OmegaCoarse, Omega);
			logger.beg("export");
				P0.exportHarwellBoeing(path + "P0.rua");
				Omega.export(path + "mesh1.ntr", { { "format", "NTR" } });
			logger.end();
		logger.end();
		logger.beg("uniform refinement");
			OmegaCoarse = Omega; // copy
			Omega.refine();
			auto P1 = FE.prolongate(OmegaCoarse, Omega);
			logger.beg("export");
				P1.exportHarwellBoeing(path + "P1.rua");
				Omega.export(path + "mesh2.ntr", { { "format", "NTR" } });
			logger.end();
		logger.end();
		logger.exp("stdin.txt");
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}