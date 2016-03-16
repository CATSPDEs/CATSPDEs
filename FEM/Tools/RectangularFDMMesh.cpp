#include "RectangularFDMMesh.hpp"

RectangularFDMMesh::RectangularFDMMesh(Point startPoint, Point endPoint, unsigned nX, unsigned nY)
	: AbstractMesh(nX * nY)
	, _horPointNumber(nX)
	, _verPointNumber(nY) 
	, _width(endPoint.x() - startPoint.x()) 
	, _height(endPoint.y() - startPoint.y())
	, _hX(_width / (nX - 1))
	, _hY(_height / (nY - 1)) {
	// TODO: throw exceptions here if needed
	for (size_t i = 0; i < nX; i++)
		for (size_t j = 0; j < nY; j++)
			_nodes[j + nX * i] = AbstractNode(startPoint.x() + i * _hX, startPoint.y() + j * _hY);
}

void RectangularFDMMesh::setBottomBC(unsigned BCType = DIRICHLET_NODE) {
	for (size_t i = 0; i < _horPointNumber; i++)
		getNode(i)->setBCType(BCType);
}

void RectangularFDMMesh::setTopBC(unsigned BCType = DIRICHLET_NODE) {
	for (size_t i = _n - 1; i >= _n - _horPointNumber; --i)
		getNode(i)->setBCType(BCType);
}

void RectangularFDMMesh::setLeftBC(unsigned BCType = DIRICHLET_NODE) {
	for (size_t i = 0; i < _n; i += _horPointNumber)
		getNode(i)->setBCType(BCType);
}

void RectangularFDMMesh::setRightBC(unsigned BCType = DIRICHLET_NODE) {
	for (size_t i = _horPointNumber - 1; i < _n; i += _horPointNumber)
		getNode(i)->setBCType(BCType);
}
