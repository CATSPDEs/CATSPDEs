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
// for global system matrix
#include "CSlCMatrix.hpp"
#include "CSCMatrix.hpp"

namespace FEM {

	namespace Mixed {

		boost::tuple<
			CSlCMatrix<double>, // diffusion + convection + reaction matrix
			CSCMatrix<double>, CSCMatrix<double>, // divirgence matrices 
			std::vector<double> // rhs vector
		>
		assembleSystem(
			OseenProblem2D const &, // (1) PDE,
			Triangulation const &, // (2) discretized domain (mesh), and BCs that connects (1) and (2):
			VectorBoundaryCondition2D const &, // natural BC,
			VectorBoundaryCondition2D const &, // essential BC;
			TriangularScalarFiniteElement const &, // for each velocity component
			TriangularScalarFiniteElement const &, // for pressure
			boost::optional<Index&> activeElementIndex = boost::none // index of the active element	
		);

	}
}
