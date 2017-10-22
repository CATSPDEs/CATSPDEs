#pragma once
#include "AbstractElement.hpp"
#include "MasterSegmentQuadratureRule.hpp"

template<Index NODE_DIM>
class SegmentElement : public AbstractElement<1, NODE_DIM, 2, 2, Node<NODE_DIM>> {
public:
	// master segment
	SegmentElement() { setX(_nodes[0], -1.); setX(_nodes[1], 1.); }
	SegmentElement(Node<NODE_DIM> const & a, Node<NODE_DIM> const & b) : AbstractElement({a, b}) {}
	AbstractQuadratureRule<1>& qRule() const override {
		return MasterSegmentQuadratureRule::instance();
	}
	virtual AbstractQuadratureRule<0>& qRuleFace() const {
		return MasterNodeQuadratureRule::instance();
	}
	std::array<Node<NODE_DIM>, 2> faces() const override {
		return { _nodes[0], _nodes[1] };
	}
	std::array<Mapping<0, 1>, 2> faceMappings() const override {
		std::array<Mapping<0, 1>, 2> M;
		M[0] = [](Void const &) { return -1.; };
		M[1] = [](Void const &) { return  1.; };
		return M;
	}
	std::array<std::array<Mapping<0, 1>, 0>, 2> faceMappingsDerivatives() const override {
		std::array<std::array<Mapping<0, 1>, 0>, 2> derM;
		return derM;
	}
	double measure() const override {
		return norm(_nodes[1] - _nodes[0]);
	}
	double diameter() const override {
		return norm(_nodes[1] - _nodes[0]);
	}
};

using SegmentElement1D = SegmentElement<1>;
using SegmentElement2D = SegmentElement<2>;
using SegmentElement3D = SegmentElement<3>;