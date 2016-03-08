#pragma once
#include"AbstractMesh.hpp"
class RectangularFDMMesh : public AbstractMesh
{
	double _width;
	double _length;
	unsigned _horDotNumber;
	unsigned _vertDotNumber;
	double _hX;
	double _hY;
public:
	unsigned nodeMapping(unsigned) const;
	RectangularFDMMesh(Point, Point, unsigned, unsigned);
};