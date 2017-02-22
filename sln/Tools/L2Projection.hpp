#pragma once
// tools
#include <boost/tuple/tuple.hpp>
// mesh
#include "Triangulation.hpp"
// matrix
#include "SymmetricCSlCMatrix.hpp" // for final linear system matrix
#include "CSCMatrix.hpp" // rhs matrix
// FEs
#include "AbstractFiniteElement.hpp"

/*
	Žilyakov Alexander, Jan 2017
*/

namespace FEM {

	boost::tuple<
		SymmetricCSlCMatrix<double>, // system mass matrix
		CSCMatrix<double> // rhs mass matrix
	>
	L2ProjectionAssembler(
		Triangulation const &, // fine and
		Triangulation const &, // coarse meshes
		TriangularScalarFiniteElement const & // FE to define interpolants
	);

}