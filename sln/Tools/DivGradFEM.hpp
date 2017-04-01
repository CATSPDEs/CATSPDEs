#pragma once
// tools
#include <boost/tuple/tuple.hpp>
// mesh
#include "Triangulation.hpp"
// PDE
#include "DiffusionReactionEqn.hpp"
#include "BoundaryCondition.hpp"
// matrix
#include "SymmetricCSlCMatrix.hpp" // for final linear system matrix
// FEs
#include "AbstractFiniteElement.hpp"

namespace FEM {

	namespace DivGrad {

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

		boost::tuple<
			SymmetricCSlCMatrix<double>, // assembled system matrix and
			std::vector<double> // rhs vector
		> assembleSystem(
			DiffusionReactionEqn2D const &, // (1) PDE,
			Triangulation const &, // (2) discretized domain (mesh), and BCs to connect (1) and (2):
			ScalarBoundaryCondition2D const &, // natural BC,
			ScalarBoundaryCondition2D const &, // essential BC
			TriangularScalarFiniteElement const &,
			boost::optional<TriangularScalarFiniteElement const &> T_FE = boost::none, // FE for (curvilinear) mapping
			boost::optional<Index&> activeElementIndex = boost::none // index of the active element	
		);

	}

}