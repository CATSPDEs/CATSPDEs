#pragma once
#include "vector.hpp"
#include "Element.hpp"
// in order to store computed images of quadrature nodes
#include "SmartMapping.hpp"

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element
// M := dimension of shapes

template <LocalIndex D, LocalIndex N, LocalIndex M>
class AbstractFiniteElement {
protected:
	double _deg;
public:
	AbstractFiniteElement(double deg) : _deg(deg) {}
	virtual ~AbstractFiniteElement() {}
	// get degree of polynomial space of shapes
	double deg() const { return _deg; };
	// get shape funcs of an element
	virtual std::vector<SmartMapping<D, M>>     getShapesOf(Element<D, N> const &) const = 0;
	// get gradients of ″
	virtual std::vector<SmartMapping<D, 2 * M>> getSGradsOf(Element<D, N> const &) const = 0;
};

	using TriangularScalarFiniteElement = AbstractFiniteElement<2, 3, 1>;
	using TriangularVectorFiniteElement = AbstractFiniteElement<2, 3, 2>;



