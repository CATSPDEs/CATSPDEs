#pragma once
#include "Point.hpp"
const unsigned INTERNAL_NODE = 0;
const unsigned DIRICHLET_NODE = 1;
const unsigned NEUMANN_NODE = 2;
const unsigned ROBIN_NODE = 3;

class AbstractNode
{
	Point _location;
	unsigned short _BC_type;
public:
	explicit AbstractNode(Point const &, unsigned BC_type = INTERNAL_NODE);
	explicit AbstractNode(double x = 0., double y = 0., double z = 0., unsigned BC_type = INTERNAL_NODE);
	void setLocation(Point const &);
	void setBC_type(unsigned);
	Point getLocation() const;
	unsigned getBC_type() const;
};