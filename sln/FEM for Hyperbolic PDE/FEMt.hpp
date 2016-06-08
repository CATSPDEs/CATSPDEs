#pragma once
// to define our problem, we need 3 objects:
// (1) PDE,
#include "HyperbolicPDE.hpp"
// (2) discretized domain (mesh),
#include "Triangulation.hpp"
// (3) and IBCs that connects (1) and (2)
#include "InitialBoundaryConditions.hpp"
#include "SymmetricContainer.hpp" // for local matrices

namespace FEMt {
	// Crank–Nicolson scheme
	vector<vector<double>> CN3(HyperbolicPDE const &, InitialBoundaryConditions const &, vector<double> const &, Triangulation&);
	// helpers
	SymmetricContainer<double> computeLocalMassMatrix(Function, array<Node, 3>&, double);
	SymmetricContainer<double> computeLocalStiffnessMatrix(Function, array<Node, 3>&, array<Node, 3>&, double);
	array<double, 3> computeLocalLoadVector(TimeFunction, double, array<Node, 3>&, array<Node, 3>&, double);
	SymmetricContainer<double> computeLocalRobinMatrix(TimeFunction, double, array<Node, 2>& nodes, double length);
	array<double, 2> computeLocalRobinVector(InitialBoundaryConditions const &, double, array<Node, 2>&, double length);
}
