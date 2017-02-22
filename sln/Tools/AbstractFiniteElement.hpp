#pragma once
// to get DOFs numn for given mesh
#include "AbstractMesh.hpp"
// in order to store computed images of quadrature nodes
#include "SmartMapping.hpp"

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element
// M := dimension of shapes

template <LocalIndex D, LocalIndex N, LocalIndex M>
class AbstractFiniteElement {
protected:
	double _deg;
	bool   _isConform;
public:
	AbstractFiniteElement(double deg, bool isConform) : _deg(deg), _isConform(isConform) {}
	virtual ~AbstractFiniteElement() {}
	// U_H := FE-space on the coarse mesh, U_h := ″ fine; return true if U_H in U_h
	bool isConform() const { return _isConform; }
	// get degree of polynomial space of shapes
	double deg() const { return _deg; };
	// get shape funcs of an element
	virtual std::vector<SmartMapping<D, M>>     getShapesOf(Element<D, N> const &) const = 0;
	// get gradients of ″
	virtual std::vector<SmartMapping<D, 2 * M>> getSGradsOf(Element<D, N> const &) const = 0;
	// DOFs
	virtual Index numbOfDOFs(AbstractMesh<D, N> const & mesh) const = 0;
	virtual std::vector<Index>   getDOFsNumeration(AbstractMesh<D, N> const & mesh, Index e) const = 0;
	virtual std::vector<Node<D>> getDOFsNodes     (AbstractMesh<D, N> const & mesh, Index e) const = 0;
	virtual std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<D, N> const & mesh, LocalIndex b) const = 0;
};

	using TriangularScalarFiniteElement = AbstractFiniteElement<2, 3, 1>;
	using TriangularVectorFiniteElement = AbstractFiniteElement<2, 3, 2>;


