#pragma once
#include"AbstractNode.hpp"
#include<vector>
#include<stdexcept>

class AbstractMesh {
protected:
	std::vector<AbstractNode> _nodes;
public:
	AbstractMesh(size_t n) : _nodes(n) {}
	virtual ~AbstractMesh() {}
	// virtual size_t nodeMapping(size_t) const = 0;
	size_t getNodesNumber() const { return _nodes.size(); };
	AbstractNode* getNode(size_t i) { return &_nodes[i]; };
};