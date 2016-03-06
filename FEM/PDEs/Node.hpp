#pragma once
#include"Point.hpp"

class Node
{
	Point _location;
	unsigned _index;
	bool _isBoundary;
public:
	Node(unsigned, Point, bool isBoundary = false);
	Node(unsigned, double, double, double, bool isBoundary = false);
	void setLocation(Point const &);
	void setIndex(unsigned);
	void setAsBoundary(bool value = true);
	Point getLocation() const;
	unsigned getIndex() const;
	bool isBoundaryNode() const;
	~Node();
};