#include "Node.hpp"

Node::Node(unsigned index, Point location, bool isBoundary) : _index(index), _location(location), _isBoundary(isBoundary)
{
}

Node::Node(unsigned index, double x, double y, double z,  bool isBoundary): _index(index), _location(x,y,z),_isBoundary(isBoundary)
{
}

void Node::setLocation(Point const & location)
{
	_location = location;
}

void Node::setIndex(unsigned index)
{
	_index = index;
}

void Node::setAsBoundary(bool value)
{
	_isBoundary = value;
}

Point Node::getLocation() const
{
	return _location;
}

unsigned Node::getIndex() const
{
	return _index;
}

bool Node::isBoundaryNode() const
{
	return _isBoundary;
}