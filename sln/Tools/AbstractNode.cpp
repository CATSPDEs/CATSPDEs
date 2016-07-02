#include "AbstractNode.hpp"

AbstractNode::AbstractNode(Point const & location, unsigned BC_type) :  _location(location), _BC_type(BC_type)
{
}

AbstractNode::AbstractNode(double x, double y, double z, unsigned BC_type): _location(x,y,z), _BC_type(BC_type)
{
}

void AbstractNode::setLocation(Point const & location)
{
	_location = location;
}

void AbstractNode::setBC_type(unsigned BC_type)
{
	_BC_type = BC_type;
}

Point AbstractNode::getLocation() const
{
	return _location;
}


unsigned AbstractNode::getBC_type() const
{
	return _BC_type;
}