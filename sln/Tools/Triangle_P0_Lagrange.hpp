#pragma once
#include "AbstractFiniteElement.hpp"
#include "Triangulation.hpp"

/*
	Alexander Žilyakov, Jan 2017
*/

class Triangle_P0_Lagrange 
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_P0_Lagrange() : AbstractFiniteElement { 0., true } {}
	Triangle_P0_Lagrange(Triangle_P0_Lagrange const &);
	Triangle_P0_Lagrange& operator=(Triangle_P0_Lagrange const &);
public:
	static auto& instance() {
		static Triangle_P0_Lagrange single;
		return single;
	}
	std::vector<SmartScalarField2D> getShapesOf(Triangle2D const &) const final {
		return { [](Node2D const &) { return 1.; } };
	}
	std::vector<SmartVectorField2D> getSGradsOf(Triangle2D const &) const final {
		return { [](Node2D const &) -> Node2D { return { 0., 0. }; } };
	}
	Index numbOfDOFs(AbstractMesh<2, 3> const & mesh) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return T->numbOfElements();
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return { e };
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Node2D> getDOFsNodes(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return { centroid(T->getElement(e)) };
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<2, 3> const & mesh, LocalIndex b) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return { 0 };
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
};