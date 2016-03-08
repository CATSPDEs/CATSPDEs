#pragma once
#include"Node.hpp"
#include<vector>
#include<stdexcept>
class AbstractMesh
{
protected:
	std::vector<Node> _nodes;
public:
	AbstractMesh();
	virtual ~AbstractMesh();
	unsigned getNodesNumber() const;
	Node getNode(unsigned) const;
	virtual unsigned nodeMapping(unsigned) const = 0;
};