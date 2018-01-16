#pragma once
#include "AbstractFiniteElement.hpp"
#include "Triangulation.hpp"

/*
	Alexander Žilyakov, Jan 2017
*/

class Triangle_P1_Lagrange 
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_P1_Lagrange() : AbstractFiniteElement { 1., true } {}
	Triangle_P1_Lagrange(Triangle_P1_Lagrange const &);
	Triangle_P1_Lagrange& operator=(Triangle_P1_Lagrange const &);
public:
	static auto& instance() {
		static Triangle_P1_Lagrange single;
		return single;
	}
	std::vector<ScalarField2D> getShapesOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const & p) { 
				return ((p[1] - t[1][1])*t[2][0] + p[0] * (t[1][1] - t[2][1]) + t[1][0] * (-p[1] + t[2][1])) / 2. / area(t); 
			},
			[=](Node2D const & p) {
				return ((-p[1] + t[0][1])*t[2][0] + t[0][0] * (p[1] - t[2][1]) + p[0] * (-t[0][1] + t[2][1])) / 2. / area(t);
			},
			[=](Node2D const & p) {
				return ((p[1] - t[0][1])*t[1][0] + p[0] * (t[0][1] - t[1][1]) + t[0][0] * (-p[1] + t[1][1])) / 2. / area(t);
			}
		};
	}
	std::vector<VectorField2D> getSGradsOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const &) { 
				return Node2D { { t[1][1] - t[2][1], -t[1][0] + t[2][0] } } / 2. / area(t);
			},
			[=](Node2D const & p) {
				return Node2D { { -t[0][1] + t[2][1], t[0][0] - t[2][0] } } / 2. / area(t);
			},
			[=](Node2D const & p) {
				return Node2D { { t[0][1] - t[1][1], -t[0][0] + t[1][0] } } / 2. / area(t);
			}
		};
	}
	Index numbOfDOFs(AbstractMesh<2, 3> const & mesh) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return T->numbOfNodes();
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodesIndicies = T->getNodesIndicies(e);
			return std::vector<Index>(nodesIndicies.begin(), nodesIndicies.end());
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Node2D> getDOFsNodes(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodes = T->getElement(e);
			return std::vector<Node2D>(nodes.begin(), nodes.end());;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<2, 3> const & mesh, LocalIndex b) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto localIndicies = excludeIndex(b);
			return std::vector<LocalIndex>(localIndicies.begin(), localIndicies.end());;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
};