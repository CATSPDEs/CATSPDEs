#include "AbstractElement.hpp"

AbstractElement::AbstractElement(unsigned index): _index(index)
{
}

AbstractElement::AbstractElement(unsigned index, std::vector<Node> const & nodes) : _nodes(nodes), _index(index)
{
}

AbstractElement::~AbstractElement()
{
}

void AbstractElement::setIndex(unsigned index)
{
	_index = index;
}

unsigned AbstractElement::getIndex()
{
	return _index;
}

Node AbstractElement::getNode(unsigned localIndex)
{
	if (localIndex < _nodes.size())
		return _nodes[localIndex];
	else
		throw std::out_of_range("local index must be in range [0,n), where n is number of nodes in element");
}

unsigned AbstractElement::getNodesNum()
{
	return _nodes.size();
}

void AbstractElement::addNode(Node const & node)
{
	_nodes.push_back(node);
}