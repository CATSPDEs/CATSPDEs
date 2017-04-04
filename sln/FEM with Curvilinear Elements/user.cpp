#include "DivGradFEM.hpp"
#include "ProjectionSolvers.hpp"

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
			Omega.refine(3);
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