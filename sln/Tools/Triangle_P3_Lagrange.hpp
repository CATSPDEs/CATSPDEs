#pragma once
#include "AbstractFiniteElement.hpp"
#include "Triangulation.hpp"

/*
	George Gubskij, Apr 2017
*/

class Triangle_P3_Lagrange
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_P3_Lagrange() : AbstractFiniteElement { 3., true } {}
	Triangle_P3_Lagrange(Triangle_P3_Lagrange const &);
	Triangle_P3_Lagrange& operator=(Triangle_P3_Lagrange const &);
public:
	static auto& instance() {
		static Triangle_P3_Lagrange single;
		return single;
	}
	std::vector<SmartScalarField2D> getShapesOf(Triangle2D const & t) const final {
		return{
			[=](Node2D const & p) {
				return 1.-5.5*p[0] + 9.*p[0] * p[0] - 4.5*p[0] * p[0] * p[0] - 5.5*p[1] + 18.*p[0]*p[1] - 13.5*p[0] * p[0] *p[1] +
					9.*p[1] * p[1] - 13.5 *p[0] *p[1] * p[1] - 4.5*p[1] * p[1] * p[1];
			},
			[=](Node2D const & p) {
				return 1. *p[0] - 4.5 *p[0]* p[0] + 4.5* p[0] * p[0] * p[0];
			},
			[=](Node2D const & p) {
				return 1.*p[1] - 4.5 *p[1] * p[1] + 4.5* p[1] * p[1] * p[1];
			},
			[=](Node2D const & p) {
				return -4.5 *p[0]* p[1] + 13.5 *p[0] * p[0]* p[1];
			},
			[=](Node2D const & p) {
				return -4.5* p[0]* p[1] + 13.5 *p[0] *p[1] * p[1];
			},
			[=](Node2D const & p) {
				return -4.5 *p[1] + 4.5* p[0]* p[1] + 18. *p[1] * p[1] - 13.5 *p[0]* p[1] * p[1] - 13.5 *p[1] * p[1] * p[1];
			},
				[=](Node2D const & p) {
				return 9.* p[1] - 22.5* p[0]* p[1] + 13.5* p[0] * p[0]* p[1] - 22.5 *p[1] * p[1] + 27. *p[0]* p[1] * p[1] + 13.5 *p[1] * p[1] * p[1];
			},
				[=](Node2D const & p) {
				return 9.* p[0] - 22.5* p[0] * p[0] + 13.5 *p[0] * p[0] * p[0] - 22.5* p[0]* p[1] + 27. *p[0] * p[0] *p[1] + 13.5* p[0] *p[1] * p[1];
			},
				[=](Node2D const & p) {
				return -4.5* p[0] + 18.* p[0] * p[0] - 13.5 *p[0] * p[0] * p[0] + 4.5* p[0]* p[1] - 13.5 *p[0] * p[0]* p[1];
			},
				[=](Node2D const & p) {
				return 27.* p[0] *p[1] - 27. *p[0] * p[0]*p[1] - 27.* p[0]* p[1] * p[1];
			}
		};
	}

	std::vector<SmartVectorField2D> getSGradsOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const & p) {
			auto x = p[0], y = p[1];
				return Node2D { { 
						-5.5 + 18. *x - 13.5 *x *x + 18.* y - 27.* x* y - 13.5 *y *y, -5.5 + 18. *x -
						13.5* x *x + 18. *y - 27.* x* y - 13.5* y*y
					} };
			},
			[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						1. - 9. *x + 13.5* x *x, 0
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						0, 1. - 9.* y + 13.5 *y *y
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						-4.5 *y + 27.* x* y, -4.5* x + 13.5* x*x
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						-4.5* y + 13.5 *y *y, -4.5* x + 27.* x *y
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						4.5* y - 13.5 *y *y, -4.5 + 4.5* x + 36. *y - 27. *x *y - 40.5 *y *y
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						-22.5* y + 27.* x *y + 27.* y *y, 9. - 22.5 *x + 13.5* x *x - 45.* y +
						54. *x *y + 40.5 *y *y
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						9. - 45. *x + 40.5* x *x - 22.5 *y + 54.* x* y + 13.5 *y *y, -22.5 *x +
						27.* x *x + 27.* x *y
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						-4.5 + 36.* x - 40.5* x *x + 4.5 *y - 27. *x *y, 4.5 *x - 13.5 *x *x
					} };
			},[=](Node2D const & p) {
				auto x = p[0], y = p[1];
				return Node2D{ {
						27. *y - 54. *x *y - 27.* y *y, 27. *x - 27. *x *x - 54.* x *y
					} };
			}
		};
	}
	Index numbOfDOFs(AbstractMesh<2, 3> const & mesh) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return T->numbOfNodes() + 2*T->numbOfRibs()+T->numbOfElements();
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodesIndicies = T->getNodesIndicies(e);
			auto ribsIndicies = T->getRibsIndicies(e);
			std::vector<Index> res(nodesIndicies.begin(), nodesIndicies.end());
			res.reserve(7);
			std::vector<Index> newRibsIndicies{
				ribsIndicies[0] * 2,ribsIndicies[0] * 2 + 1,
				ribsIndicies[1] * 2,ribsIndicies[1] * 2 + 1,
				ribsIndicies[2] * 2,ribsIndicies[2] * 2 + 1,
				e + 2 * T->numbOfRibs()
			};
			for (auto& el : newRibsIndicies) el += T->numbOfNodes();
			res.insert(res.end(), newRibsIndicies.begin(), newRibsIndicies.end());
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Node2D> getDOFsNodes(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto enodes = T->getElement(e);
			std::vector<Node2D> mnodes;
			mnodes.reserve(6);
			for (LocalIndex i : {0, 1, 2})
			{
				mnodes.emplace_back(T->getRibNode(e, i, 1. / 3.));
				mnodes.emplace_back(T->getRibNode(e, i, 2. / 3.));
			}
			std::vector<Node2D> res;
			res.reserve(10);
			res.insert(res.end(), enodes.begin(), enodes.end());
			res.insert(res.end(), mnodes.begin(), mnodes.end());
			res.emplace_back(centroid(enodes));
			return res;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<2, 3> const & mesh, LocalIndex b) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			if (b == 0) return{ 1,3,4,2 };
			if (b == 1) return{ 2,5,6,0};
			return{ 0,7,8,1};
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
};