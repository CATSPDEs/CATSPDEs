#pragma once
#include"Node.hpp"
#include<vector>
#include<stdexcept>

const unsigned INDEX_IS_NOT_USED = 0;

class AbstractElement
{
protected:
	std::vector<Node> _nodes;
	unsigned _index;
public:
	AbstractElement(unsigned index=INDEX_IS_NOT_USED);
	AbstractElement(unsigned, std::vector<Node> const &);
	virtual ~AbstractElement();
	void setIndex(unsigned);
	void addNode(Node const &);
	unsigned getIndex();
	Node getNode(unsigned localIndex); //localIndex [0,n) is the local index of node in element
	unsigned getNodesNumber();
};