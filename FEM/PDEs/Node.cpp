#include "Node.hpp"

Node::Node(unsigned index, Point location, bool isBoundary) : _index(index), _location(location), _isBoundary(isBoundary)
{
}

Node::Node(unsigned index, REAL x, REAL y, REAL z,  bool isBoundary): _index(index), _location(x,y,z),_isBoundary(isBoundary)
{
}

void Node::setLocation(Point location)
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