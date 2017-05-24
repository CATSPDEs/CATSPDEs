#pragma once
#include "AbstractMesh.hpp"

/*
	Alexander Žilyakov, March 2017
*/

class SegmentMesh
	: public AbstractMesh<1, 2> {
	// Inherited via AbstractMesh
public:
	// redefine
	Index numbOfElements() const final { return _nodes.size() - 1; }
	std::array<Index, 2> getNodesIndicies(Index e) const { return { e, e + 1 }; }
	Segment1D getElement(Index e) const final { return { _nodes[e], _nodes[e + 1] }; }
	// pure virtual
	SegmentMesh& refine(Index numbOfRefinements = 1) final;
	std::vector<Index> getFineNeighborsIndicies(Index e) const final;
	SegmentMesh& import(std::istream & from = std::cin) final;
	void export(std::ostream & to = std::cout, Parameters const & param = {}) const final;
	// in order to work w/ strings, not streams
	using AbstractMesh::import;
	using AbstractMesh::export;
};