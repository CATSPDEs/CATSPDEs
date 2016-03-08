#include "AbstractMesh.hpp"

AbstractNode AbstractMesh::getNode(unsigned localIndex) const
{
	unsigned pureLocalIndex = nodeMapping(localIndex);
	return _nodes[pureLocalIndex];
}
