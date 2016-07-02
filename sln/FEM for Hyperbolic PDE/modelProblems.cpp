#include <fstream>
#include <string>
#include "FEM_t.hpp"

// input R × R × T —> R functions for IBCs
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
	try {
		HyperbolicPDE PDE(chi, sigma, a, f);
		InitialConditions ICs(initialPosition, initialVelocity);
		// BCs
		DirichletBC_t Dirichlet(G1_D, G1);
		NeumannBC_t Neumann(G2_N, G2);
		RobinBC_t Robin(G3_R, G3_N, G3);
		BoundaryConditions_t BCs({ &Dirichlet, &Neumann, &Robin });
		TimeFrames time(0., 1., 1); // 0.___.5___.1 — initially we have 3 time frames ([0., 1.] refined once)
		vector<vector<double>> soln;
		string path("Mathematica/Model Problem Analysis/tests/"), currentPath;
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
			vector<vector<double>> uVec(soln.size(), soln[0]);
			for (size_t m = 0; m < time.size(); ++m)
				for (size_t j = 0; j < Omega.numbOfNodes(); ++j)
					uVec[m][j] = u(Omega.getNode(j), time[m]);
			oU << uVec;
			// save meshes parameters
			vector<double> edges = Omega.longestEdges();
			oH << *max_element(edges.begin(), edges.end());
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