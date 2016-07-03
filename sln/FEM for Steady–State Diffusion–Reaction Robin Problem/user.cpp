#include <fstream>
#include <string>
#include "FEM.hpp"

// input R × R —> R functions (PDE and BCs) and the mesh (domain)
//#include "PureRobinProblem.hpp"
//#include "PureNeumannProblem.hpp"
//#include "PureDirichletProblem.hpp"
#include "kittyProblem.hpp"

int main() {
	try {
		DiffusionReactionEqn PDE(a, c, f);
		BoundaryConditions BCs(ListOfBCs);
		vector<double> soln(Omega.numbOfNodes());
		string path("Mathematica/Model Problem Analysis/");
		ofstream output;
		unsigned refCount;
		cout << "numb of refinements: ";
		cin >> refCount;
		for (unsigned i = 0; i <= refCount; ++i) {
			soln = FEM::computeDiscreteSolution(PDE, Omega, BCs);
			// save the result
			output.open(path + "xi" + to_string(i) + ".dat");
			output << soln;
			output.close();
			Omega.save(ofstream(path + "n" + to_string(i) + ".dat"), ofstream(path + "t" + to_string(i) + ".dat"));
			Omega.refine();
		}
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}