#pragma once
#include "Triangle_P1_Lagrange.hpp"

/*
	Alexander Žilyakov, Apr 2017
*/

class Triangle_Pt3_LagrangeBubble 
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_Pt3_LagrangeBubble() : AbstractFiniteElement { 3., true } {}
	Triangle_Pt3_LagrangeBubble(Triangle_Pt3_LagrangeBubble const &);
	Triangle_Pt3_LagrangeBubble& operator=(Triangle_Pt3_LagrangeBubble const &);
public:
	static auto& instance() {
		static Triangle_Pt3_LagrangeBubble single;
		return single;
	}
	std::vector<SmartScalarField2D> getShapesOf(Triangle2D const & t) const final {
		auto shapes = Triangle_P1_Lagrange::instance().getShapesOf(t);
		auto bubble = [=](Node2D const & p) {
			return 27. * shapes[0](p) * shapes[1](p) * shapes[2](p);
		};
		shapes.emplace_back(bubble);
		return shapes;
	}
	std::vector<SmartVectorField2D> getSGradsOf(Triangle2D const & t) const final {
		auto shapes = Triangle_P1_Lagrange::instance().getShapesOf(t);
		auto sGrads = Triangle_P1_Lagrange::instance().getSGradsOf(t);
		auto bGrad = [=](Node2D const & p) {
			return 27. * (
				shapes[1](p) * shapes[2](p) * sGrads[0](p) +
				shapes[0](p) * shapes[2](p) * sGrads[1](p) +
				shapes[0](p) * shapes[1](p) * sGrads[2](p)
			);
		};
		sGrads.emplace_back(bGrad);
		return sGrads;
	}
	Index numbOfDOFs(AbstractMesh<2, 3> const & mesh) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return T->numbOfNodes() + T->numbOfElements();
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodesIndicies = T->getNodesIndicies(e);
			std::vector<Index> res;
			res.reserve(4);
			res.insert(res.end(), nodesIndicies.begin(), nodesIndicies.end());
			res.emplace_back(T->numbOfNodes() + e);
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Node2D> getDOFsNodes(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodes = T->getElement(e);
			std::vector<Node2D> res;
			res.reserve(4);
			res.insert(res.end(), nodes.begin(), nodes.end());
			res.emplace_back(centroid(nodes));
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<2, 3> const & mesh, LocalIndex b) const final {
		return Triangle_P1_Lagrange::instance().getBndryDOFsLocalIndicies(mesh, b);
	}
};