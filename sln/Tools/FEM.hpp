#pragma once
// to define our problem, we need 3 objects:
// (1) PDE,
#include "DiffusionReactionEqn.hpp"
// (2) discretized domain (mesh),
#include "Triangulation.hpp"
// (3) and BCs that connects (1) and (2)
#include "BoundaryConditions.hpp"
#include "SymmetricContainer.hpp" // for local matrices

// our model problem:
//
// –nabla . (a nabla u) + cu = f,             if (x, y) in Omega,
// –n . (a nabla u) = kappa (u – g_D) – g_N,  if (x, y) in bndry of Omega,
//
// where a > 0, c >= c_0 > 0, kappa > 0, f, g_D, and g_N are given R × R —> R functions
//
// ( we will get input R × R —> R functions a, c, and f from PDE data structure (1), 
//   kappa, g_D, and g_N—from BCs data strucrure (2), and
//   Omega—from mesh data strucrure (3) )
//
// we convert our problem into a discrete one:
//
// (massMatrix + stiffnessMatrix + robinMatrix) . xi = loadVector + robinVector 
// (or A.xi = b for short),
//
// where matrix / element entries are given:
//
// massMatrix(i, j)      = dintt_Omega { c hatFunction_i hatFunction_j },
// stiffnessMatrix(i, j) = dintt_Omega { a nabla hatFunction_i . nabla hatFunction_j },
// robinMatrix(i, j)     = dintt_{bndry of Omega} { kappa hatFunction_i hatFunction_j },
// loadVector(i)         = dintt_Omega { f hatFunction_i },
// robinVector(i)        = dintt_{bndry of Omega} { (kappa g_D + g_N) hatFunction_i },
//
// where hatFunction_i denotes linear basis function taking unity on ith node and zero elsewhere
// so we have to assemble and solve n × n linear system, n := numb of nodes of the mesh Omega

namespace FEM {
	// classic finite element method
	vector<double> computeDiscreteSolution(DiffusionReactionEqn const &, Triangulation&, BoundaryConditions&);
	// helpers
	SymmetricContainer<double> computeLocalMassMatrix     (Function, array<Node, 3>&, double);
	SymmetricContainer<double> computeLocalStiffnessMatrix(Function, array<Node, 3>&, array<Node, 3>&, double);
	array<double, 3>           computeLocalLoadVector     (Function, array<Node, 3>&, array<Node, 3>&, double);
	SymmetricContainer<double> computeLocalRobinMatrix    (BoundaryConditions const &, array<Node, 2>&, double);
	array<double, 2>           computeLocalRobinVector    (BoundaryConditions const &, array<Node, 2>&, double);
}
