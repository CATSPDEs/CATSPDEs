#pragma once
#include "AbstractQuadratureRule.hpp"

template<
	LocalIndex ELEM_DIM, 
	LocalIndex NODE_DIM, 
	LocalIndex NUMB_OF_NODES,
	LocalIndex NUMB_OF_FACES,
	typename FACE_TYPE
>
class AbstractElement {
protected:
	std::array<Node<NODE_DIM>, NUMB_OF_NODES> _nodes;
public:
	// todo: https://stackoverflow.com/questions/301203/extract-c-template-parameters
	typedef FACE_TYPE faceType;
	AbstractElement() { 
		static_assert(NODE_DIM >= ELEM_DIM, "invalid node / elem dimension");
		for (auto& p : _nodes) zeroOut(p); 
	}
	AbstractElement(std::array<Node<NODE_DIM>, NUMB_OF_NODES> const & iniList) : _nodes(iniList) {
		static_assert(NODE_DIM >= ELEM_DIM, "invalid node / elem dimension");
	}
	virtual ~AbstractElement() {}
	auto  nodes() const { return _nodes; }
	auto& nodes() { return _nodes; }
	auto numbOfFaces() const { return NUMB_OF_FACES; }
	virtual AbstractQuadratureRule<ELEM_DIM>& qRule() const = 0;
	virtual AbstractQuadratureRule<ELEM_DIM - 1>& qRuleFace() const {
		FACE_TYPE faceElem;
		return faceElem.qRule();
	}
	virtual std::array<FACE_TYPE, NUMB_OF_FACES> faces() const = 0;
	virtual std::array<Mapping<ELEM_DIM - 1, ELEM_DIM>, NUMB_OF_FACES> faceMappings() const = 0;
	virtual std::array<std::array<Mapping<ELEM_DIM - 1, ELEM_DIM>, ELEM_DIM - 1>, NUMB_OF_FACES> faceMappingsDerivatives() const = 0;
	virtual double measure() const = 0;
	virtual double diameter() const = 0;
};
