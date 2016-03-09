#pragma once
#include"AbstractNode.hpp"
#include<vector>
#include<stdexcept>

class AbstractMesh
{
protected:
	size_t _n; // numb of nodes
	std::vector<AbstractNode> _nodes;
public:
	AbstractMesh(size_t n) : _n(n), _nodes(n) {}
	virtual ~AbstractMesh() {}
	unsigned getNodesNumber() const { return _n; };
	AbstractNode getNode(unsigned) const;
	virtual unsigned nodeMapping(unsigned) const = 0;
};