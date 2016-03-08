#include "AbstractNode.hpp"

AbstractNode::AbstractNode(Point location, unsigned BCType) :  _location(location), _BCType(BCType)
{
}

AbstractNode::AbstractNode(double x, double y, double z, unsigned BCType): _location(x,y,z), _BCType(BCType)
{
}

void AbstractNode::setLocation(Point const & location)
{
	_location = location;
}

void AbstractNode::setBCType(unsigned BCType)
{
	_BCType = BCType;
}

Point AbstractNode::getLocation() const
{
	return _location;
}


unsigned AbstractNode::getBCType() const
{
	return _BCType;
}