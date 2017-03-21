#pragma once
#include "AbstractFiniteElement.hpp"
#include "SegmentMesh.hpp"

/*
	Alexander Žilyakov, March 2017
*/

class Segment_P1_Lagrange 
	: public SegmentFiniteElement {
	// singleton
	Segment_P1_Lagrange() : AbstractFiniteElement { 1., true } {}
	Segment_P1_Lagrange(Segment_P1_Lagrange const &);
	Segment_P1_Lagrange& operator=(Segment_P1_Lagrange const &);
public:
	static auto& instance() {
		static Segment_P1_Lagrange single;
		return single;
	}
	std::vector<SmartScalarField1D> getShapesOf(Segment1D const & s) const final {
		return {
			[=](double const & p) { return (s[1] - p) / length(s); },
			[=](double const & p) { return (p - s[0]) / length(s); }
		};
	}
	std::vector<SmartScalarField1D> getSGradsOf(Segment1D const & s) const final {
		throw std::logic_error("not implemented");
	}
	Index numbOfDOFs(AbstractMesh<1, 2> const & mesh) const final {
		return mesh.numbOfNodes();
	}
	std::vector<Index> getDOFsNumeration(AbstractMesh<1, 2> const & mesh, Index e) const final {
		return { e, e + 1 };
	}
	std::vector<double> getDOFsNodes(AbstractMesh<1, 2> const & mesh, Index e) const final {
		return { mesh.getNode(e), mesh.getNode(e + 1) };
	}
	std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<1, 2> const & mesh, LocalIndex b) const final {
		throw std::logic_error("not implemented");
	}
};