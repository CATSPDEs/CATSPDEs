#include <fstream>
#include <string>
#include "FEMt.hpp"

// input R × R —> R and R × R × T —> R functions (PDE and IBCs)
#include "allBCsExample.hpp"
// and the mesh (domain)
#include "Triangulation.hpp"

int main() {
	try {
		HyperbolicPDE PDE(chi, sigma, a, f);
		InitialBoundaryConditions IBCs(g_D, g_N, kappa, initialVelocity);
		Triangulation Omega(generateMesh());
		Omega.refine(2);
		vector<double> time = { 1., 2., 3., 4., 5., 6., 7. };
		vector<vector<double>> soln = FEMt::CN3(PDE, IBCs, time, Omega);

		string path("Mathematica/Model Problem Analysis/");
		ofstream output;
		
		Omega.save(ofstream(path + "n.dat"), ofstream(path + "t.dat"));
		output.open(path + "xi.dat");
		output << soln;

		//vector<double> soln(Omega.numbOfNodes());
		//string path("Mathematica/Model Problem Analysis/");
		//ofstream output;
		//unsigned refCount;
		//cout << "numb of refinements: ";
		//cin >> refCount;
		//for (unsigned i = 0; i < refCount; ++i) {
		//	Omega.refine();
		//	soln = FEM::computeDiscreteSolution(PDE, BCs, Omega);
		//	// save the result
		//	output.open(path + "xi" + to_string(i) + ".dat");
		//	output << soln;
		//	output.close();
		//	Omega.save(ofstream(path + "n" + to_string(i) + ".dat"), ofstream(path + "t" + to_string(i) + ".dat"));
		//}
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}