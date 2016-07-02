#pragma once
// to define our problem, we need 3 objects:
// (1) PDE,
#include "HyperbolicPDE.hpp"
// (2a) discretized domain (mesh),
#include "Triangulation.hpp"
// (2b) discretized time interval (time frames),
#include "TimeFrames.hpp"
// (3) and IBCs that connects (1) and (2)
#include "InitialConditions.hpp"
#include "BoundaryConditions_t.hpp"
// for local matrices
#include "SymmetricContainer.hpp" 

namespace FEM_t {
	// Crank–Nicolson scheme, 3 time frames
	vector<vector<double>> CN3(HyperbolicPDE const &, TimeFrames const &, Triangulation&, InitialConditions const &, BoundaryConditions_t&, 
		                       Function_t u = emptyFunc_t); // exact soln (if we are dealing w / model problem)
	// Backward difference formula, 3 time frames
	vector<vector<double>> BDF3(HyperbolicPDE const &, TimeFrames const &, Triangulation&, InitialConditions const &, BoundaryConditions_t&, 
		                        Function_t u = emptyFunc_t);
	// helpers
	SymmetricContainer<double> computeLocalMassMatrix     (Function, array<Node, 3>&, double);
	SymmetricContainer<double> computeLocalStiffnessMatrix(Function, array<Node, 3>&, array<Node, 3>&, double);
	array<double, 3>           computeLocalLoadVector     (Function_t, double, array<Node, 3>&, array<Node, 3>&, double);
	SymmetricContainer<double> computeLocalRobinMatrix    (BoundaryConditions_t const &, double, array<Node, 2>&, double length);
	array<double, 2>           computeLocalRobinVector    (BoundaryConditions_t const &, double, array<Node, 2>&, double length);
}
