#include <fstream>
#include <string>
#include "FEM.hpp"

inline bool G1(Node& p) { // left half of kitty’s face
	if (p.x() <= 0.) return true;
	return false;
}

inline double G1_D(Node& p) {
	return p.x() * p.x() + p.y() * p.y();
}

int main() {
	ifstream iNodes("Mathematica/Generate Mesh/nodes.dat"),
	         iTriangles("Mathematica/Generate Mesh/triangles.dat");
	ofstream oXi("Mathematica/xi.dat"),
	         oDirichletNodes("Mathematica/DirichletNodes.dat"),
	         oNodes("Mathematica/nodes.dat"), 
	         oTriangles("Mathematica/triangles.dat");
	try {
		// import kitty mesh
		Triangulation Omega(iNodes, iTriangles);
		Omega.refine();
		//DiffusionReactionEqn LaplaceEqn;
		//BoundaryConditions BCs({
		//	make_shared<DirichletBC>(G1_D, G1),
		//	make_shared<NeumannBC>() // free bndry on right half of kitty’s face
		//});
		//// solve
		//vector<double> soln = FEM::computeDiscreteSolution(LaplaceEqn, Omega, BCs);
		//// save the result
		//oXi << soln;
		Boundary bndry = Omega.computeBoundary();
		oDirichletNodes << FEM::computeBoundaryNodes(Omega, bndry, G1);
		Omega.save(oNodes, oTriangles);
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}