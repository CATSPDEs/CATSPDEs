#include "SingletonLogger.hpp"
// FEs to build prolongation matrices for
#include "Triangle_P0_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
#include "Triangle_P2_Lagrange.hpp"
// prolongation matrix data structure
#include "CSCMatrix.hpp"

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
				&Triangle_P2_Lagrange::instance()
			};
			auto FEIndex = logger.opt("choose finite element", { "Lagrange P0", "Lagrange P1", "Lagrange P2" });
			auto& FE = *FEs[FEIndex];
		logger.end();
		logger.beg("import mesh");
			Triangulation Omega;
			Omega.import(path + "crs_uniform.ntn");
			Omega.enumerateRibs();
		logger.end();
		logger.beg("non-uniform refinement");
			Index n;
			logger.inp("numb of elements to refine", n);
			vector<Index> elems2refine_untreated(n);
			logger.inp("indicies of elements to refine", elems2refine_untreated);
			Indicies elems2refine(elems2refine_untreated.begin(), elems2refine_untreated.end());
			auto OmegaCoarse = Omega; // copy
			Omega.refine(elems2refine);
			CSCMatrix<double> P0(FE.numbOfDOFs(OmegaCoarse), FE.numbOfDOFs(Omega));
			P0.generatePatternFrom(createDOFsConnectivityList(OmegaCoarse, Omega, FE));
			FE.prolongate(P0, OmegaCoarse, Omega);
			logger.beg("export");
				P0.exportHarwellBoeing(path + "P0.rua");
				Omega.export(path + "mesh1.ntr", { { "format", "NTR" } });
			logger.end();
		logger.end();

		logger.exp("stdin.txt");
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}