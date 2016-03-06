#pragma once
#include"Node.hpp"
#include<vector>
#include<stdexcept>
class AbstractMesh
{
protected:
	std::vector<Node> _nodes;
	std::vector<Node> _boundaryNodes;
public:
	AbstractMesh();
	virtual ~AbstractMesh();
	unsigned getNodesNumber() const;
	unsigned getBoundaryNodesNumber() const;
	Node getNode(unsigned) const;
	virtual unsigned nodeMapping(unsigned) const = 0;
};