#include "AbstractMesh.hpp"

AbstractMesh::AbstractMesh()
{
}

AbstractMesh::~AbstractMesh()
{
}

unsigned AbstractMesh::getNodesNumber() const
{
	return _nodes.size();
}

unsigned AbstractMesh::getBoundaryNodesNumber() const
{
	return _boundaryNodes.size();
}

Node AbstractMesh::getNode(unsigned localIndex) const
{
	unsigned pureLocalIndex = nodeMapping(localIndex);
	return _nodes[pureLocalIndex];
}

