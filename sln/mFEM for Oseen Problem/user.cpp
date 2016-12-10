#include "SingletonLogger.hpp"
#include "FEM.hpp"

// logger
SingletonLogger& logger = SingletonLogger::instance();

using namespace FEM::Mixed;
using std::string;

int main() {
	string iPath("Mathematica/");
	try {
		// set PDE
		OseenProblem2D PDE(
			1.,
			[](Node2D const & p) -> Node2D {
				return	{ p[0], p[1] };
			},
			.5,
			[](Node2D const & p) -> Node2D {
				return { 1., 1. };
			},
			[](Node2D const & p) {
				return norm(p);
			}
		);
		// import and set up mesh
		Triangulation Omega;
		Omega.AbstractMesh::import(iPath + "mesh.ntn");
		Omega.enumerateRibs();
		Omega.AbstractMesh::export(iPath + "mesh.ntr", { {"format", "NTR"} });
		// set up BCs
		DirichletVectorCondition2D DirichletBC(
			[](Node2D const &) -> Node2D {
				return { 0., 0. };
			}
		);

		P2P1Assembler(PDE, Omega, DirichletBC);

		//logger.buf << PDE.massTerm();
		//logger.log();
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}