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
	size_t nodeMapping(size_t i) const { return i; } // here we have simple mapping
	RectangularFDMMesh(Point, Point, unsigned, unsigned);
	void setBottomBC(unsigned);
	void setTopBC(unsigned);
	void setLeftBC(unsigned);
	void setRightBC(unsigned);
};