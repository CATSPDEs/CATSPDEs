#include "SegmentMesh.hpp"

SegmentMesh& SegmentMesh::refine(Index numbOfRefinements) {
	throw std::logic_error("not implemented");
}

std::vector<Index> SegmentMesh::getFineNeighborsIndicies(Index e) const {
	throw std::logic_error("not implemented");
}

SegmentMesh& SegmentMesh::import(std::istream & from) {
	Index n;
	from >> n;
	_nodes.resize(n);
	from >> _nodes;
	return *this;
}

void SegmentMesh::export(std::ostream & to, Parameters const & param) const {
	throw std::logic_error("not implemented");
}
