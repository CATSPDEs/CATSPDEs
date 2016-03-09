#pragma once
#include"AbstractNode.hpp"
#include<vector>
#include<stdexcept>

const unsigned INDEX_IS_NOT_USED = 0;

class AbstractElement
{
protected:
	std::vector<AbstractNode> _nodes;
	unsigned _index;
public:
	AbstractElement(unsigned index=INDEX_IS_NOT_USED);
	AbstractElement(unsigned, std::vector<AbstractNode> const &);
	virtual ~AbstractElement();
	void setIndex(unsigned);
	void addNode(AbstractNode const &);
	unsigned getIndex();
	AbstractNode getNode(unsigned localIndex); //localIndex [0,n) is the local index of node in element
	unsigned getNodesNumber();
};