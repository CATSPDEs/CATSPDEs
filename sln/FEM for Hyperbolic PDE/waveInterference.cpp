#include <fstream>
#include <string>
#include "FEM_t.hpp"
#include "Triangulation.hpp"

// input R × R × T —> R functions for BCs

inline bool G1(Node& p) {
	if (p.x() == -3.) return true;
	return false;
}
double const PI = 3.14159265358979323846;
inline double G1_D(Node& p, double t) {
	return .0125 * sin(1.5 * PI * t); // sine wave!!
}

inline bool G2(Node& p) {
	if (p.x() * p.x() + p.y() * p.y() <= 1.21) return true;
	return false;
}

int main() {
	unsigned spaceRefCount, timeRefCount;
	cout << "refine SPACE [...] times: ";
	cin >> spaceRefCount;
	cout << "refine TIME [...] times: ";
	cin >> timeRefCount;
	try {
		HyperbolicPDE waveEqn; // hom. wave eqn	
		InitialConditions ICs; // zero position and velocity
		// BCs
		DirichletBC_t DirichletWave(G1_D, G1); // sine wave at strips’ ends,
		DirichletBC_t DirichletHom (zeroFunc_t, G2); // fixed hole, and
		NeumannBC_t Neumann; // free boundary elsewhere  
		BoundaryConditions_t BCs({ &DirichletWave, &DirichletHom, &Neumann });
		TimeFrames time(0., 15., timeRefCount); // time frames
		Triangulation Omega(generateMesh()); // mesh w/ “strips” and hole—to demonstrate wave inerference
		Omega.refine(spaceRefCount);
		// sovle w/ Crank–Nicolson scheme
		vector<vector<double>> soln = FEM_t::CN3(waveEqn, time, Omega, ICs, BCs);
		// save the result
		string path("Mathematica/Wave Interference Animation/");
		ofstream oTime(path + "time.dat"), 
			     oXi(path + "xi.dat"), 
			     oG1(path + "G1.dat"),
			     oG2(path + "G2.dat");
		time.save(oTime);
		oXi << soln;
		Boundary bndry = Omega.computeBoundary();
		for (size_t j = 0; j < bndry.size(); ++j) {
			if (G1(Omega.getNode(bndry[j][0]))) oG1 << bndry[j][0] << '\n';
			else if (G2(Omega.getNode(bndry[j][0]))) oG2 << bndry[j][0] << '\n';
		}
		Omega.save(ofstream(path + "nodes.dat"), ofstream(path + "triangles.dat"));
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}