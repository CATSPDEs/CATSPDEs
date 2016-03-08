#pragma once
#include"AbstractNode.hpp"
#include<vector>
#include<stdexcept>
class AbstractMesh
{
protected:
	std::vector<AbstractNode> _nodes;
public:
	AbstractMesh();
	virtual ~AbstractMesh();
	unsigned getNodesNumber() const;
	AbstractNode getNode(unsigned) const;
	virtual unsigned nodeMapping(unsigned) const = 0;
};