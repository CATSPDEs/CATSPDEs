#include <fstream>
#include <string>
#include "FEM.hpp"

// model soln
inline double u(Node& p) {
	return p.x() * p.x() + p.y() * p.y() + 1.;
}

// input R × R —> R functions for PDE
inline double a(Node& p) {
	return 1.;
}
inline double c(Node& p) {
	return 3.;
}
inline double f(Node& p) {
	return  3. * ( p.x() * p.x() + p.y() * p.y() ) - 1.;
}

// BCs

// pure Neumann
inline double N1(Node& p) {
	return 2.;
}
// pure Robin
inline double N2(Node& p) {
	return 4.;
}

int main() {
	char problem = '\0';
	while (problem != 'D' && problem != 'N' && problem != 'R') {
		cout << "choose pure BCs (D, N or R): ";
		cin >> problem;
	}
	unsigned refCount;
	cout << "numb of refinements: ";
	cin >> refCount;
	string path("Mathematica/Model Problem Analysis/tests/"), currentPath;
	try {
		Triangulation Omega; // unit square mesh
		DiffusionReactionEqn PDE(a, c, f);
		// BCs
		shared_ptr<AbstractBC> BC;
		if      (problem == 'D') BC = make_shared<DirichletBC>(u);       // pure Dirichlet test,
		else if (problem == 'N') BC = make_shared<NeumannBC>(N1);        // pure Neumann, and
		else                     BC = make_shared<RobinBC>(oneFunc, N2); // pure Robin
		BoundaryConditions BCs({ BC });
		vector<double> soln(Omega.numbOfNodes());
		for (unsigned i = 0; i < refCount; ++i) {
			// (1) refine
			Omega.refine();
			// (2) solve
			soln = FEM::computeDiscreteSolution(PDE, Omega, BCs);
			// (3) save the result
			currentPath = path + to_string(i + 1) + "/";
			ofstream oXi(currentPath + "xi.dat"),
			         oU(currentPath + "u.dat"),
			         oH(currentPath + "deltaH.dat");
			oXi << soln;
			oU << FEM::constructVector(u, Omega); // save vector constructed from model soln
			oH << Omega.longestEdge();
			Omega.save(ofstream(currentPath + "nodes.dat"), ofstream(currentPath + "triangles.dat"));
		}
		ofstream oNumbOfTests(path + "numbOfRows.dat");
		oNumbOfTests << refCount;
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}