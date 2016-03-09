#pragma once
#include"AbstractMesh.hpp"

class RectangularFDMMesh : public AbstractMesh
{
	double _width;
	double _height;
	unsigned _horPointNumber;
	unsigned _verPointNumber;
	double _hX;
	double _hY;
public:
	unsigned nodeMapping(unsigned) const;
	RectangularFDMMesh(Point, Point, unsigned, unsigned);
};