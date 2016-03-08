#include "Node.hpp"

Node::Node(unsigned index, Point location, unsigned BCType) : _index(index), _location(location), _BCType(BCType)
{
}

Node::Node(unsigned index, double x, double y, double z, unsigned BCType): _index(index), _location(x,y,z), _BCType(BCType)
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

void Node::setBCType(unsigned BCType)
{
	_BCType = BCType;
}

Point Node::getLocation() const
{
	return _location;
}

unsigned Node::getIndex() const
{
	return _index;
}

unsigned Node::getBCType() const
{
	return _BCType;
}