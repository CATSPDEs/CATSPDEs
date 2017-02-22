#include <fstream>
#include "FEM.hpp"
#include "SingletonLogger.hpp"

// model soln
inline double u(Node const & p) {
	return p.x() * p.x() + p.y() * p.y() + 1.;
}

// input R × R —> R functions for PDE
inline double a(Node const & p) {
	return 1.;
}
inline double c(Node const & p) {
	return 3.;
}
inline double f(Node const & p) {
	return  3. * ( p.x() * p.x() + p.y() * p.y() ) - 1.;
}

// BCs

// pure Neumann
inline double N1(Node const & p) {
	return 2.;
}
// pure Robin
inline double N2(Node const & p) {
	return 4.;
}



int main() {
	string oPath("Mathematica/Model Problem Analysis/tests/"), oCurrentPath;
	// logger
	SingletonLogger& logger = SingletonLogger::instance();
	try {
		// choose problem (pure BCs)
		size_t problemIndex = logger.opt(
			"choose pure BCs",
			{ "Dirichlet", "Neumann", "Robin" }
		); 

		// numb of mesh refinements
		logger.inp("numb of refinements");
		unsigned refCount;
		cin >> refCount;

		// set the problem
		logger.beg("generate coarse mesh, set PDE, set BCs");
			Triangulation Omega; // unit square mesh
			DiffusionReactionEqn PDE(a, c, f);
			// BCs
			shared_ptr<AbstractBC> BC;
			if (problemIndex == 0)      BC = make_shared<DirichletBC>(u);       // pure Dirichlet test,
			else if (problemIndex == 1) BC = make_shared<NeumannBC>(N1);        // pure Neumann, and
			else                        BC = make_shared<RobinBC>(oneFunc, N2); // pure Robin
			BoundaryConditions BCs({ BC });
			vector<double> soln(Omega.numbOfNodes());
		logger.end();

		// solve
		for (unsigned i = 0; i < refCount; ++i) {
			logger.beg("iteration " + to_string(i + 1));

				// (1) refine
				logger.beg("refine mesh");
					Omega.refine();
				logger.end();

				// (2) solve
				logger.beg("compute FE-solution");
					soln = FEM::computeDiscreteSolution(PDE, Omega, BCs);
				logger.end();

				// (3) save the result
				logger.beg("save results to " + oPath);
					oCurrentPath = oPath + to_string(i + 1) + "/";
					ofstream oXi(oCurrentPath + "xi.dat"),
					         oH (oCurrentPath + "h.dat");
					oXi << soln;
					oH << Omega.longestEdge();
					Omega.save(ofstream(oCurrentPath + "nodes.dat"), ofstream(oCurrentPath + "triangles.dat"));
				logger.end();

			logger.end();
		}

		ofstream oRefCount    (oPath + "refCount.dat"),
		         oProblemIndex(oPath + "problemIndex.dat");
		oRefCount << refCount;
		oProblemIndex << problemIndex;
	}
	catch (exception const & e) {
		logger.err(e.what());
	}
	return 0;
}