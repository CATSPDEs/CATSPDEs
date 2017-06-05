#pragma once
#include "AbstractFiniteElement.hpp"
#include "Triangulation.hpp"
// prolongation
#include "CSCMatrix.hpp"

/*
	Alexander Žilyakov, Jan 2017
*/

class Triangle_P1_CrouzeixRaviart
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_P1_CrouzeixRaviart() : AbstractFiniteElement { 1., false } {}
	Triangle_P1_CrouzeixRaviart(Triangle_P1_CrouzeixRaviart const &);
	Triangle_P1_CrouzeixRaviart& operator=(Triangle_P1_CrouzeixRaviart const &);
public:
	static auto& instance() {
		static Triangle_P1_CrouzeixRaviart single;
		return single;
	}
	std::vector<SmartScalarField2D> getShapesOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const & p) {
				return ((-2 * p[1] + t[0][1] + t[1][1])*t[2][0] + t[1][0] * (2 * p[1] - t[0][1] - t[2][1]) - (2 * p[0] - t[0][0])*(t[1][1] - t[2][1])) / 2. / area(t);
			},
			[=](Node2D const & p) {
				return ((2 * p[1] - t[0][1] - t[1][1])*t[2][0] + (2 * p[0] - t[1][0])*(t[0][1] - t[2][1]) + t[0][0] * (-2 * p[1] + t[1][1] + t[2][1])) / 2. / area(t);
			},
			[=](Node2D const & p) {
				return 1. + (2 * ((-p[1] + t[0][1])*t[1][0] + t[0][0] * (p[1] - t[1][1]) + p[0] * (-t[0][1] + t[1][1]))) / 2. / area(t);
			}
		};
	}
	std::vector<SmartVectorField2D> getSGradsOf(Triangle2D const & t) const final {
		return {
			[=](Node2D const &) { 
				return Node2D { { t[2][1] - t[1][1], t[1][0] - t[2][0] } } / area(t);
			},
			[=](Node2D const & p) {
				return Node2D { { t[0][1] - t[2][1], -t[0][0] + t[2][0] } } / area(t);
			},
			[=](Node2D const & p) {
				return Node2D { { t[1][1] - t[0][1], t[0][0] - t[1][0] } } / area(t);
			}
		};
	}
	Index numbOfDOFs(AbstractMesh<2, 3> const & mesh) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return T->numbOfRibs();
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto ribsIndicies = T->getRibsIndicies(e);
			return std::vector<Index>(ribsIndicies.begin(), ribsIndicies.end());
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<Node2D> getDOFsNodes(AbstractMesh<2, 3> const & mesh, Index e) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) {
			auto nodes = midNodes(T->getElement(e));
			return std::vector<Node2D>(nodes.begin(), nodes.end());;
		}
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<2, 3> const & mesh, LocalIndex b) const final {
		auto T = dynamic_cast<Triangulation const *>(&mesh);
		if (T) return { b };
		throw std::invalid_argument("FE interpolant is not defined on this mesh");
	}
	void prolongate(CSCMatrix<double>& P, AbstractMesh<2, 3> const & cMesh, AbstractMesh<2, 3> const & fMesh) const final {
		std::unordered_map<Index, Index> freq;
		auto chop = [](double x) { return fabs(x) < 1e-10 ? 0. : x; };
		for (Index ci = 0; ci < cMesh.numbOfElements(); ++ci) {
			auto shapes = getShapesOf(cMesh.getElement(ci));
			auto nodes = getDOFsNodes(fMesh, ci);
			auto rows = getDOFsNumeration(cMesh, ci),
			     cols = getDOFsNumeration(fMesh, ci);
			for (LocalIndex j = 0; j < cols.size(); ++j) {
				++freq[cols[j]];
				for (LocalIndex i = 0; i < rows.size(); ++i)
					P(rows[i], cols[j]) += .5 * chop(shapes[i](nodes[j]));
			}
			for (Index fi : fMesh.getFineNeighborsIndicies(ci)) {
				cols = getDOFsNumeration(fMesh, fi);
				nodes = getDOFsNodes(fMesh, fi);
				for (LocalIndex j = 0; j < cols.size(); ++j) {
					++freq[cols[j]];
					for (LocalIndex i = 0; i < rows.size(); ++i)
						P(rows[i], cols[j]) += .5 * chop(shapes[i](nodes[j]));
				}
			}
		}			
		for (auto const & kvp : freq) if (kvp.second == 1) P.modifyColumn(kvp.first, [&](double& val) { return val *= 2.; });
	}
};