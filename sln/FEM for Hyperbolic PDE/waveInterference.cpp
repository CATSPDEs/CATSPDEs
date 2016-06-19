#include <fstream>
#include <string>
#include "FEMt.hpp"
#include "Triangulation.hpp"

#include "array.hpp"

// input R × R × T —> R functions for BCs

inline double RobinCoefficient(Node& p, double t) {
	if (p.x() == -3. ||
	    p.x() * p.x() + p.y() * p.y() <= 1.21) return 10e50;
	return 0.;
}

double const PI = 3.14159265358979323846;

inline double DirichletCondition(Node& p, double t) {
	if (p.x() == -3.) return .0125 * sin(1.5 * PI * t); // sine wave!!
	return 0.;
}

int main() {
	try {
		HyperbolicPDE waveEqn; // hom. wave eqn		
		InitialBoundaryConditions IBCs(DirichletCondition, zeroTimeFunc, RobinCoefficient); // zero initial velocity by default
		Triangulation Omega(generateMesh()); // mesh w/ “strips” and hole—to demonstrate wave inerference
		Omega.refine(2);
		vector<double> time; // time frames
		for (double t = 0.; t <= 10.; t += .025)
			time.push_back(t);
		// sovle w/ Crank–Nicolson scheme
		vector<vector<double>> soln = FEMt::CN3(waveEqn, IBCs, time, Omega);
		// save the result
		string path("Mathematica/Model Problem Analysis/");
		ofstream oTime(path + "time.dat"), oXi(path + "xi.dat"), oBndry(path + "bndry.dat");
		oTime << time;
		oXi << soln;
		oBndry << Omega.computeBoundary();
		Omega.save(ofstream(path + "n.dat"), ofstream(path + "t.dat"));
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}