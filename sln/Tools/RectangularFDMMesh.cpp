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

void RectangularFDMMesh::setBottomBC(unsigned BC_type = DIRICHLET_NODE) {
	for (size_t i = 0; i < _horPointNumber; i++)
		getNode(i)->setBC_type(BC_type);
}

void RectangularFDMMesh::setTopBC(unsigned BC_type = DIRICHLET_NODE) {
	for (size_t i = _nodes.size() - 1; i >= _nodes.size() - _horPointNumber; --i)
		getNode(i)->setBC_type(BC_type);
}

void RectangularFDMMesh::setLeftBC(unsigned BC_type = DIRICHLET_NODE) {
	for (size_t i = 0; i < _nodes.size(); i += _horPointNumber)
		getNode(i)->setBC_type(BC_type);
}

void RectangularFDMMesh::setRightBC(unsigned BC_type = DIRICHLET_NODE) {
	for (size_t i = _horPointNumber - 1; i < _nodes.size(); i += _horPointNumber)
		getNode(i)->setBC_type(BC_type);
}
