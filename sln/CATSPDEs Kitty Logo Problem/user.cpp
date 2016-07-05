#include <fstream>
#include <string>
#include "FEM.hpp"

inline bool fixedBoundary(Node& p) { // left half of kitty’s face
	if (p.x() <= 0.) return true;
	return false;
}

inline double f(Node& p) {
	return p.norm();
}

int main() {
	ifstream iNodes("Mathematica/Generate Mesh/nodes.dat"),
	         iTriangles("Mathematica/Generate Mesh/triangles.dat");
	ofstream oXi("Mathematica/Draw Logo/xi.dat"),
	         oDirichletNodes("Mathematica/Draw Logo/DirichletNodes.dat"),
	         oNodes("Mathematica/Draw Logo/nodes.dat"), 
	         oTriangles("Mathematica/Draw Logo/triangles.dat");
	try {
		// import kitty mesh
		Triangulation Omega(iNodes, iTriangles);
		DiffusionReactionEqn PoissonEqn(oneFunc, zeroFunc, f);
		BoundaryConditions BCs({
			make_shared<DirichletBC>(oneFunc, fixedBoundary),
			make_shared<NeumannBC>() // free bndry on right half of kitty’s face
		});
		// solve
		vector<double> soln = FEM::computeDiscreteSolution(PoissonEqn, Omega, BCs);
		// save the result
		oXi << soln;
		Boundary bndry = Omega.computeBoundary();
		oDirichletNodes << FEM::computeBoundaryNodes(Omega, bndry, fixedBoundary);
		Omega.save(oNodes, oTriangles);
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}