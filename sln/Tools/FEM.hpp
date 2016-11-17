#pragma once
#include <boost/tuple/tuple.hpp>
#include "CSCMatrix.hpp" // restriction (prolongation) operators
#include "SymmetricCSlCMatrix.hpp" // for final linear system matrix

#include "DiffusionReactionEqn.hpp"
#include "Triangulation.hpp"
#include "BoundaryConditions.hpp"

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

	namespace DivGrad {

		boost::tuple<
			SymmetricCSlCMatrix<double>, // assembled system matrix and
			std::vector<double> // rhs vector
		> linearLagrangeAssembler(
			DiffusionReactionEqn2D const &, // (1) PDE,
			Triangulation const &, // (2) discretized domain (mesh),
			DirichletScalarCondition2D const & // (3) and BCs that connects (1) and (2)
		);

		namespace Multigrid {

			extern std::function<std::vector<double>(
				SymmetricCSlCMatrix<double>&, // system matrix
				std::vector<double> const &, // rhs
				std::vector<double> const & // initial guess
			)> smoother;
			extern Index gamma; // numb of recursive coarse iterations (1 for V–cycle, 2 for W–cycle)

			boost::tuple<
				SymmetricCSlCMatrix<double>, // assembled system matrix and
				std::vector<double> // rhs vector
			> linearLagrangeSetter(
				DiffusionReactionEqn2D const & PDE,
				DirichletScalarCondition2D const & DirichletBC,
				Triangulation&, // initial coarse mesh
				Index // numb of mesh levels
			);

			std::vector<double> iteration( // for solving A.z = f
				Index, // current mesh level
				std::vector<double> const &, // initial approximation to soln
				std::vector<double> const & // rhs
			);

		}

	}
	
	namespace Mixed {
		
		//boost::tuple<
		//	SymmetricCSlCMatrix<double>, // assembled system matrix and
		//	std::vector<double> // rhs vector
		//> linearLagrangeAssembler(
		//	OseenProblem2D const &, // (1) PDE,
		//	Triangulation const &, // (2) discretized domain (mesh),
		//	DirichletScalarCondition2D const & // (3) and BCs that connects (1) and (2)
		//);

	}

	//// local matrices and vectors
	//SymmetricContainer<double> computeLocalMassMatrix     (Function, array<Node, 3>&, double);
	//SymmetricContainer<double> computeLocalStiffnessMatrix(Function, array<Node, 3>&, array<Node, 3>&, double);
	//array<double, 3>           computeLocalLoadVector     (Function, array<Node, 3>&, array<Node, 3>&, double);
	//SymmetricContainer<double> computeLocalRobinMatrix    (BoundaryConditions const &, array<Node, 2>&, double);
	//array<double, 2>           computeLocalRobinVector    (BoundaryConditions const &, array<Node, 2>&, double);
	//// make vector of model soln
	//vector<double> constructVector(Function, Triangulation&);
	//// save boundary nodes
	//vector<size_t> computeBoundaryNodes(Triangulation&, Boundary&, Predicate);
}
