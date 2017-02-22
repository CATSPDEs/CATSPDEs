#include "L2Projection.hpp"
#include "ProjectionSolvers.hpp" // conjugate gradients

// finite elements to test L2 projection with
#include "Triangle_P0_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_CrouzeixRaviart.hpp"

using std::vector;
using std::string;
using std::to_string;
using namespace FEM;
using namespace ProjectionSolvers::Krylov;

int main() {
	string path("Mathematica/");
	auto& logger = SingletonLogger::instance();
	try {
		if (logger.opt("what shall I do", { "prepare meshes", "compute interpolants" })) {
			logger.beg("import meshes and fine dofs");
				Triangulation cMesh, fMesh;
				cMesh.import(path + "cMesh.ntnr");
				fMesh.import(path + "fMesh.ntnr");
			logger.end();
			auto& FE = Triangle_P1_CrouzeixRaviart::instance();
			logger.beg("import fine dofs");
				vector<double> fDOFs(FE.numbOfDOFs(fMesh));
				import(fDOFs, path + "fDOFs.dat");
			logger.end();
			logger.beg("run L2 assembler");
				auto matrices = L2ProjectionAssembler(fMesh, cMesh, FE);
				auto& lhs = boost::get<0>(matrices);
				auto& rhs = boost::get<1>(matrices);
			logger.end();
			logger.beg("compute L2 projection w/ CG");
				auto cDOFs = boost::get<0>(CG(lhs, rhs * fDOFs));
			logger.end();
			logger.beg("export matrices");
				static_cast<CSCMatrix<double>>(static_cast<CSlCMatrix<double>>(lhs)).exportHarwellBoeing(path + "lhs.rsa");
				rhs.exportHarwellBoeing(path + "rhs.rra");
			logger.end();
			logger.beg("export coarse dofs");
				export(cDOFs, path + "cDOFs.dat");
			logger.end();
		}
		else {
			Triangulation Omega;
			logger.beg("import NTN-mesh");
				Omega.import(path + "mesh.ntn");
			logger.end();
			logger.beg("enum ribs, export coarse and fine NTNR meshes");
				Omega.enumerateRibs().export(path + "cMesh.ntnr", { {"format", "NTNR"} });
				Omega.refine().export(path + "fMesh.ntnr", { { "format", "NTNR" } });
			logger.end();
		}
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}