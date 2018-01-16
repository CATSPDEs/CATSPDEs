#pragma once
#include "AbstractFiniteElement.hpp"
#include "Triangulation.hpp"

/*
	Alexander Žilyakov, Jan 2017
*/

class Triangle_P2_Lagrange
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_P2_Lagrange() : AbstractFiniteElement { 2., true } {}
	Triangle_P2_Lagrange(Triangle_P2_Lagrange const &);
	Triangle_P2_Lagrange& operator=(Triangle_P2_Lagrange const &);
public:
	static auto& instance() {
		static Triangle_P2_Lagrange single;
		return single;
	}
	std::vector<ScalarField2D> getShapesOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const & p) {
				return (((-2.*p[1] + t[0][1] + t[1][1])*t[2][0] + t[1][0] * (2.*p[1] - t[0][1] - t[2][1]) - (2.*p[0] - t[0][0])*(t[1][1] - t[2][1]))*
					   ((-p[1] + t[1][1])*t[2][0] + t[1][0] * (p[1] - t[2][1]) + p[0] * (-t[1][1] + t[2][1]))) / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return (((-2.*p[1] + t[0][1] + t[1][1])*t[2][0] - (2.*p[0] - t[1][0])*(t[0][1] - t[2][1]) + t[0][0] * (2.*p[1] - t[1][1] - t[2][1]))*
					   ((-p[1] + t[0][1])*t[2][0] + t[0][0] * (p[1] - t[2][1]) + p[0] * (-t[0][1] + t[2][1]))) / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return (((-p[1] + t[0][1])*t[1][0] + t[0][0] * (p[1] - t[1][1]) + p[0] * (-t[0][1] + t[1][1]))*(-(t[0][1] - t[1][1])*(2.*p[0] - t[2][0]) +
					   t[0][0] * (2.*p[1] - t[1][1] - t[2][1]) + t[1][0] * (-2.*p[1] + t[0][1] + t[2][1]))) / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return (4.*((p[1] - t[0][1])*t[1][0] + p[0] * (t[0][1] - t[1][1]) + t[0][0] * (-p[1] + t[1][1]))*
				       ((-p[1] + t[0][1])*t[2][0] + t[0][0] * (p[1] - t[2][1]) + p[0] * (-t[0][1] + t[2][1]))) / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return (4.*((-p[1] + t[0][1])*t[1][0] + t[0][0] * (p[1] - t[1][1]) + p[0] * (-t[0][1] + t[1][1]))*
				       ((-p[1] + t[1][1])*t[2][0] + t[1][0] * (p[1] - t[2][1]) + p[0] * (-t[1][1] + t[2][1]))) / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return (4.*((p[1] - t[1][1])*t[2][0] + p[0] * (t[1][1] - t[2][1]) + t[1][0] * (-p[1] + t[2][1]))*
				       ((-p[1] + t[0][1])*t[2][0] + t[0][0] * (p[1] - t[2][1]) + p[0] * (-t[0][1] + t[2][1]))) / 4. / pow(area(t), 2);
			}
		};
	}
	std::vector<VectorField2D> getSGradsOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const & p) {
				return Node2D { { 
						(t[1][1] - t[2][1])*(t[1][1] * (4.*p[0] - t[0][0] - 3.*t[2][0]) + t[0][1] * (t[1][0] - t[2][0]) + 4.*p[1] * (-t[1][0] + t[2][0]) + (-4.*p[0] + t[0][0] + 3.*t[1][0])*t[2][1]),
						(t[1][0] - t[2][0])*(4.*p[1] * (t[1][0] - t[2][0]) + t[0][1] * (-t[1][0] + t[2][0]) + t[1][1] * (-4.*p[0] + t[0][0] + 3.*t[2][0]) + (4.*p[0] - t[0][0] - 3.*t[1][0])*t[2][1])
					} } / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return Node2D { {
						(t[0][1] - t[2][1])*(-t[0][1] * t[1][0] - (3.*t[0][1] + t[1][1])*t[2][0] + 4.*p[1] * (-t[0][0] + t[2][0]) + 4.*p[0] * (t[0][1] - t[2][1]) + t[1][0] * t[2][1] + t[0][0] * (t[1][1] + 3.*t[2][1])),
						(t[0][0] - t[2][0])*(t[0][1] * (-4.*p[0] + t[1][0]) + 4.*p[1] * (t[0][0] - t[2][0]) + (3.*t[0][1] + t[1][1])*t[2][0] + 4.*p[0] * t[2][1] - t[1][0] * t[2][1] - t[0][0] * (t[1][1] + 3.*t[2][1]))
					} } / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return Node2D { {
						(t[0][1] - t[1][1])*(4.*p[1] * (-t[0][0] + t[1][0]) + 4.*p[0] * (t[0][1] - t[1][1]) + t[1][1] * (3.*t[0][0] + t[2][0]) - t[0][1] * (3.*t[1][0] + t[2][0]) + (t[0][0] - t[1][0])*t[2][1]),
						(t[0][0] - t[1][0])*(4.*p[1] * (t[0][0] - t[1][0]) + 4.*p[0] * (-t[0][1] + t[1][1]) - t[1][1] * (3.*t[0][0] + t[2][0]) + t[0][1] * (3.*t[1][0] + t[2][0]) + (-t[0][0] + t[1][0])*t[2][1])
					} } / 4. / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return Node2D { {
						(t[0][1] * (p[0] - t[1][0]) + p[1] * (-t[0][0] + t[1][0]) + (-p[0] + t[0][0])*t[1][1])*(-t[0][1] + t[2][1]) + (t[0][1] - t[1][1])*(p[1] * (t[0][0] - t[2][0]) + t[0][1] * (-p[0] + t[2][0]) + (p[0] - t[0][0])*t[2][1]),
						(t[0][1] * (p[0] - t[1][0]) + p[1] * (-t[0][0] + t[1][0]) + (-p[0] + t[0][0])*t[1][1])*(t[0][0] - t[2][0]) + (-t[0][0] + t[1][0])*(p[1] * (t[0][0] - t[2][0]) + t[0][1] * (-p[0] + t[2][0]) + (p[0] - t[0][0])*t[2][1])
					} } / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return Node2D { {
						(p[1] * (t[0][0] - t[1][0]) + t[0][1] * (-p[0] + t[1][0]) + (p[0] - t[0][0])*t[1][1])*(-t[1][1] + t[2][1]) + (-t[0][1] + t[1][1])*(p[1] * (t[1][0] - t[2][0]) + t[1][1] * (-p[0] + t[2][0]) + (p[0] - t[1][0])*t[2][1]),
						(p[1] * (t[0][0] - t[1][0]) + t[0][1] * (-p[0] + t[1][0]) + (p[0] - t[0][0])*t[1][1])*(t[1][0] - t[2][0]) + (t[0][0] - t[1][0])*(p[1] * (t[1][0] - t[2][0]) + t[1][1] * (-p[0] + t[2][0]) + (p[0] - t[1][0])*t[2][1])
					} } / pow(area(t), 2);
			},
			[=](Node2D const & p) {
				return Node2D { { 
						(t[1][1] - t[2][1])*(p[1] * (t[0][0] - t[2][0]) + t[0][1] * (-p[0] + t[2][0]) + (p[0] - t[0][0])*t[2][1]) + (-t[0][1] + t[2][1])*(t[1][1] * (p[0] - t[2][0]) + p[1] * (-t[1][0] + t[2][0]) + (-p[0] + t[1][0])*t[2][1]), 
						(-t[1][0] + t[2][0])*(p[1] * (t[0][0] - t[2][0]) + t[0][1] * (-p[0] + t[2][0]) + (p[0] - t[0][0])*t[2][1]) + (t[0][0] - t[2][0])*(t[1][1] * (p[0] - t[2][0]) + p[1] * (-t[1][0] + t[2][0]) + (-p[0] + t[1][0])*t[2][1])
					} } / pow(area(t), 2);
			}
		};
	}
	Index numbOfDOFs(AbstractMesh<2, 3> const & mesh) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return T->numbOfNodes() + T->numbOfRibs();
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodesIndicies = T->getNodesIndicies(e);
			auto ribsIndicies = T->getRibsIndicies(e);
			auto n = T->numbOfNodes();
			for (auto& el : ribsIndicies) el += n;
			std::vector<Index> res(nodesIndicies.begin(), nodesIndicies.end());
			res.reserve(3);
			res.insert(res.end(), ribsIndicies.begin(), ribsIndicies.end());
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Node2D> getDOFsNodes(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto enodes = T->getElement(e);
			auto mnodes = midNodes(enodes);
			std::vector<Node2D> res;
			res.reserve(6);
			res.insert(res.end(), enodes.begin(), enodes.end());
			res.insert(res.end(), mnodes.begin(), mnodes.end());
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<2, 3> const & mesh, LocalIndex b) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto localIndicies = excludeIndex(b);
			std::vector<LocalIndex> res;
			res.reserve(3);
			res.insert(res.end(), localIndicies.begin(), localIndicies.end());
			res.emplace_back(3 + b);
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
};