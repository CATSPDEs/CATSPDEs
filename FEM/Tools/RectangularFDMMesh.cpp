#include "RectangularFDMMesh.hpp"

unsigned RectangularFDMMesh::nodeMapping(unsigned) const
{
	return 0;  //to be declared
}

RectangularFDMMesh::RectangularFDMMesh(Point startPoint, Point endPoint, unsigned nX, unsigned nY)
	: AbstractMesh(nX * nY)
	, _horPointNumber(nX)
	, _verPointNumber(nY) 
	, _width(abs(endPoint.x() - startPoint.x())) 
	, _height(abs(endPoint.y() - startPoint.y()))
	, _hX(_width / nX)
	, _hY(_height / nY) {
	// TODO: throw exceptions here if needed
	for (size_t i = 0; i < nX; i++)
		for (size_t j = 0; j < nY; j++)
			_nodes[j + nX * i] = AbstractNode(startPoint.x() + i * _hX, startPoint.y() + j * _hY);
}
