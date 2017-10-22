#pragma once
#include "SegmentElement.hpp"
#include "MasterTriangleQuadratureRule.hpp"

template<Index NODE_DIM>
class TriangleElement : public AbstractElement<2, NODE_DIM, 3, 3, SegmentElement<NODE_DIM>> {
public:
	// master triangle
	TriangleElement() {
		_nodes[1][0] = _nodes[2][1] = 1.;
	}
	TriangleElement(Node<NODE_DIM> const & a, Node<NODE_DIM> const & b, Node<NODE_DIM> const & c) : AbstractElement({ a, b, c }) {}
	AbstractQuadratureRule<2>& qRule() const override {
		return MasterTriangleQuadratureRule::instance();
	}
	std::array<SegmentElement<NODE_DIM>, 3> faces() const override {
		return { 
			SegmentElement<NODE_DIM>(_nodes[1], _nodes[2]),
			SegmentElement<NODE_DIM>(_nodes[2], _nodes[0]),
			SegmentElement<NODE_DIM>(_nodes[0], _nodes[1])
		};
	}
	std::array<Mapping<1, 2>, 3> faceMappings() const override {
		std::array<Mapping<1, 2>, 3> M;
		M[0] = [](double const & p) -> Node2D { return { .5 * (1. - p), .5 * (1. + p) }; };
		M[1] = [](double const & p) -> Node2D { return { 0., .5 * (1. - p) }; };
		M[2] = [](double const & p) -> Node2D { return { .5 * (1. + p), 0. }; };
		return M;
	}
	std::array<std::array<Mapping<1, 2>, 1>, 3> faceMappingsDerivatives() const override {
		std::array<std::array<Mapping<1, 2>, 1>, 3> derM;
		derM[0][0] = [](double const & p) -> Node2D { return { -.5, .5 }; };
		derM[1][0] = [](double const & p) -> Node2D { return { 0., -.5 }; };
		derM[2][0] = [](double const & p) -> Node2D { return { .5, 0. };  };
		return derM;
	}
	double measure() const override {
		return norm(crossProduct(_nodes[1] - _nodes[0], _nodes[2] - _nodes[0])) / 2.;
	}
	double diameter() const override { // longest edge
		auto edges = faces();
		std::array<double, 3> lengths;
		std::transform(edges.begin(), edges.end(), lengths.begin(), [](auto const & edge) { return edge.measure(); });
		return *std::max_element(lengths.begin(), lengths.end());
	}
};

using TriangleElement2D = TriangleElement<2>;
using TriangleElement3D = TriangleElement<3>;