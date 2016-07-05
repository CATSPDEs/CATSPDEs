#include <fstream>
#include <string>
#include "FEM_t.hpp"

// input R × R × T —> R functions for PDE and IBCs
#include "sineWave.hpp"
//#include "simple.hpp"

int main() {
	unsigned numbOfTests, spaceRefCount, timeRefCount;
	string method;
	cout << "numb of tests (rows in computeConvergenceAnalysisTables.nb): ";
	cin >> numbOfTests;
	cout << "each time refine SPACE [...] times: ";
	cin >> spaceRefCount;
	cout << "each time refine TIME [...] times: ";
	cin >> timeRefCount;
	while (method != "CN3" && method != "BDF3") {
		cout << "method (CN3, BDF3): ";
		cin >> method;
	}
	string path("Mathematica/Model Problem Analysis/tests/"), currentPath;
	try {
		Triangulation Omega(Node(0., 0.), Node(1., 1.)); // simple square mesh
		HyperbolicPDE PDE(chi, sigma, a, f);
		InitialConditions ICs(initialPosition, initialVelocity);
		BoundaryConditions_t BCs({
			make_shared<DirichletBC_t>(G1_D, G1),
			make_shared<NeumannBC_t>(G2_N, G2),
			make_shared<RobinBC_t>(G3_R, G3_N)
		});
		TimeFrames time(0., 1., 1); // 0.___.5___.1 — initially we have 3 time frames ([0., 1.] refined once)
		vector<vector<double>> soln;
		for (unsigned i = 0; i < numbOfTests; ++i) {
			// (1) solve
			soln = method == "CN3" ? FEM_t::CN3 (PDE, time, Omega, ICs, BCs, u) 
			                       : FEM_t::BDF3(PDE, time, Omega, ICs, BCs, u);
			// (2) save the results
			currentPath = path + to_string(i + 1) + "/";
			ofstream oTime(currentPath + "time.dat"),
				     oXi(currentPath + "xi.dat"),
				     oU(currentPath + "u.dat"),
				     oH(currentPath + "deltaH.dat"),
				     oT(currentPath + "deltaT.dat"),
				     oDirichletNodes(currentPath + "DirichletNodes.dat"),
				     oNeumannNodes  (currentPath + "NeumannNodes.dat"),
				     oRobinNodes    (currentPath + "RobinNodes.dat");
			time.save(oTime);
			oXi << soln;
			Omega.save(ofstream(currentPath + "nodes.dat"), ofstream(currentPath + "triangles.dat"));
			// save model soln
			oU << FEM_t::constructVector(u, time, Omega);
			// save meshes parameters
			vector<double> edges = Omega.longestEdges();
			oH << Omega.longestEdge();
			oT << time[1] - time[0]; // const delta t
			// for analysis
			Boundary bndry = Omega.computeBoundary();
			for (size_t j = 0; j < bndry.size(); ++j) {
				if (G1(Omega.getNode(bndry[j][0]))) oDirichletNodes << bndry[j][0] << '\n';
				else if (G2(Omega.getNode(bndry[j][0]))) oNeumannNodes << bndry[j][0] << '\n';
				else oRobinNodes << bndry[j][0] << '\n';
			}
			// refine
			Omega.refine(spaceRefCount);
			time.refine(timeRefCount);
		}	
		ofstream oNumbOfTests(path + "numbOfRows.dat");
		oNumbOfTests << numbOfTests;
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}