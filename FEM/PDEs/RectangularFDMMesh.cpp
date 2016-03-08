#include "RectangularFDMMesh.hpp"

unsigned RectangularFDMMesh::nodeMapping(unsigned) const
{
	return 0;  //to be declared
}

RectangularFDMMesh::RectangularFDMMesh(Point startPoint, Point endPoint, unsigned nX, unsigned nY):_horDotNumber(nX), _vertDotNumber(nY)
{
	_width = abs(endPoint.x - startPoint.x);
	_length = abs(endPoint.y - startPoint.y);
	_hY = _length / nY;
	_hX = _width / nX;
	_nodes.reserve(nX*nY);
	for (int i = 0; i < nX; i++)
	{
		for (int j; j < nY; j++)
		{
			AbstractNode node(startPoint.x + i*_hX, startPoint.y + j*_hY, 0);
			_nodes.push_back(node);
		}
	}
}
