#pragma once
#include "AbstractFiniteElement.hpp"
#include "SegmentMesh.hpp"

/*
	Alexander Žilyakov, June 2017
*/

class Segment_P3_Hermite 
	: public SegmentFiniteElement {
	// singleton
	Segment_P3_Hermite() : AbstractFiniteElement { 3., true } {}
	Segment_P3_Hermite(Segment_P3_Hermite const &);
	Segment_P3_Hermite& operator=(Segment_P3_Hermite const &);
public:
	static auto& instance() {
		static Segment_P3_Hermite single;
		return single;
	}
	std::vector<SmartScalarField1D> getShapesOf(Segment1D const & s) const final {
		auto denom = pow(s[0], 4.) - 4 * pow(s[0], 3.) * s[1] + 6. * pow(s[0], 2.) * pow(s[1], 2.) - 4. * s[0] * pow(s[1], 3.) + pow(s[1], 4.);
		return {
			[=](double const & p) { 
				return (
					(-2. * s[0] + 2. * s[1]) * pow(p, 3.) +
					(3. * pow(s[0], 2.) - 3. * pow(s[1], 2.)) * pow(p, 2.) +
					(-6. * pow(s[0], 2.) * s[1] + 6. * s[0] * pow(s[1], 2.)) * p +
					3. * pow(s[0], 2.) * pow(s[1], 2.) - 4. * s[0] * pow(s[1], 3.) + pow(s[1], 4.)
				) / denom;
			},
			[=](double const & p) {
				return (
					(pow(s[0], 2) - 2 * s[0]*s[1] + pow(s[1], 2)) * pow(p, 3.) +
					(-pow(s[0], 3) + 3 * s[0]*pow(s[1], 2) - 2 * pow(s[1], 3)) * pow(p, 2.) +
					(2 * pow(s[0], 3)*s[1] - 3 * pow(s[0], 2)*pow(s[1], 2) + pow(s[1], 4)) * p -
					(pow(s[0], 3)*pow(s[1], 2)) + 2 * pow(s[0], 2)*pow(s[1], 3) - s[0]*pow(s[1], 4)
				) / denom;
			},
			[=](double const & p) {
				return (
					(2 * s[0] - 2 * s[1]) * pow(p, 3.) +
					(-3 * pow(s[0], 2) + 3 * pow(s[1], 2)) * pow(p, 2.) +
					(6 * pow(s[0], 2)*s[1] - 6 * s[0]*pow(s[1], 2)) * p +
					pow(s[0], 4) - 4 * pow(s[0], 3)*s[1] + 3 * pow(s[0], 2)*pow(s[1], 2)
				) / denom;
			},
			[=](double const & p) {
				return (
					(pow(s[0], 2) - 2 * s[0]*s[1] + pow(s[1], 2)) * pow(p, 3.) +
					(-2 * pow(s[0], 3) + 3 * pow(s[0], 2)*s[1] - pow(s[1], 3)) * pow(p, 2.) +
					(pow(s[0], 4) - 3 * pow(s[0], 2)*pow(s[1], 2) + 2 * s[0]*pow(s[1], 3)) * p -
					(pow(s[0], 4)*s[1]) + 2 * pow(s[0], 3)*pow(s[1], 2) - pow(s[0], 2)*pow(s[1], 3)
				) / denom;
			}
		};
	}
	std::vector<SmartScalarField1D> getSGradsOf(Segment1D const & s) const final {
		throw std::logic_error("not implemented");
	}
	Index numbOfDOFs(AbstractMesh<1, 2> const & mesh) const final {
		return 2 * mesh.numbOfNodes();
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<1, 2> const & mesh, Index e) const final {
		return { e, e + mesh.numbOfNodes(), e + 1, e + 1 + mesh.numbOfNodes() };
	}
	std::vector<double> getDOFsNodes(AbstractMesh<1, 2> const & mesh, Index e) const final {
		return { mesh.getNode(e), mesh.getNode(e), mesh.getNode(e + 1), mesh.getNode(e + 1) };
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<1, 2> const & mesh, LocalIndex b) const final {
		throw std::logic_error("not implemented");
	}
};