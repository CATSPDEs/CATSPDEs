#pragma once
// tools
#include <boost/tuple/tuple.hpp>
#pragma once
// mesh
#include "Triangulation.hpp"
// PDE
#include "OseenProblem.hpp"
#include "BoundaryCondition.hpp"
// FEs
#include "AbstractFiniteElement.hpp"

namespace FEM {

	namespace Mixed {

		//boost::tuple<
		//	SymmetricCSlCMatrix<double>, // assembled system matrix and
		//	std::vector<double> // rhs vector
		//> 
		void assembleSystem(
			OseenProblem2D const &, // (1) PDE,
			Triangulation const &, // (2) discretized domain (mesh), and BCs that connects (1) and (2):
			VectorBoundaryCondition2D const &, // essential BC, 
			VectorBoundaryCondition2D const &, // natural BC;
			TriangularScalarFiniteElement const &, // for each velocity component
			TriangularScalarFiniteElement const & // for pressure
			);

	}
}
