#pragma once
#include"Point.hpp"
const unsigned INTERNAL_NODE = 0;
const unsigned DIRICHLET_NODE = 1;
const unsigned NEUMANN_NODE = 2;
const unsigned ROBIN_NODE = 3;

class Node
{
	Point _location;
	unsigned _index;
	unsigned short _BCType;
public:
	Node(unsigned, Point, unsigned BCType=INTERNAL_NODE);
	Node(unsigned, double, double, double, unsigned BCType = INTERNAL_NODE);
	void setLocation(Point const &);
	void setIndex(unsigned);
	void setBCType(unsigned);
	Point getLocation() const;
	unsigned getIndex() const;
	unsigned getBCType() const;
};