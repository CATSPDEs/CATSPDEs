#pragma once
#include "SegmentElement.hpp"

template<Index NODE_DIM>
class QuadrilateralElement : public AbstractElement<2, NODE_DIM, 4, 4, SegmentElement<NODE_DIM>> {
public:
	// master square
	QuadrilateralElement() {
		_nodes[2][0] = _nodes[2][1] = _nodes[1][0] = _nodes[3][1] =  1.;
		_nodes[0][0] = _nodes[0][1] = _nodes[1][1] = _nodes[3][0] = -1.;
	}
	QuadrilateralElement(Node<NODE_DIM> const & a, Node<NODE_DIM> const & b, Node<NODE_DIM> const & c, Node<NODE_DIM> const & d) : AbstractElement({ a, b, c, d }) {}
	AbstractQuadratureRule<2>& qRule() const override {
		throw std::logic_error("not implemented");
	}
	std::array<SegmentElement<NODE_DIM>, 4> faces() const override {
		return {
			SegmentElement<NODE_DIM>(_nodes[0], _nodes[1]),
			SegmentElement<NODE_DIM>(_nodes[1], _nodes[2]),
			SegmentElement<NODE_DIM>(_nodes[2], _nodes[3]),
			SegmentElement<NODE_DIM>(_nodes[3], _nodes[0])
		};
	}
	std::array<Mapping<1, 2>, 4> faceMappings() const override {
		std::array<Mapping<1, 2>, 4> M;
		M[0] = [](double const & p) -> Node2D { return { .5 * (1. - p), .5 * (1. + p) }; };
		M[1] = [](double const & p) -> Node2D { return { 0., .5 * (1. - p) }; };
		M[2] = [](double const & p) -> Node2D { return { .5 * (1. + p), 0. }; };
		return M;
	}
	std::array<std::array<Mapping<1, 2>, 1>, 4> faceMappingsDerivatives() const override {
		std::array<std::array<Mapping<1, 2>, 1>, 4> M;
		M[0][0] = [](double const & p) -> Node2D { return { -.5, .5 }; };
		M[1][0] = [](double const & p) -> Node2D { return { 0., -.5 }; };
		M[2][0] = [](double const & p) -> Node2D { return { .5, 0. };  };
		return M;
	}
	double measure() const override {
		return .5 * (
			norm(crossProduct(_nodes[1] - _nodes[0], _nodes[2] - _nodes[0])) +
			norm(crossProduct(_nodes[2] - _nodes[0], _nodes[3] - _nodes[0]))
		);
	}
	double diameter() const override { // longest diag
		std::vector<SegmentElement<NODE_DIM>> edges { 
			{ _nodes[0], _nodes[2] }, 
			{ _nodes[1], _nodes[3] } 
		};
		std::array<double, 2> lengths;
		std::transform(edges.begin(), edges.end(), lengths.begin(), [](auto const & edge) { return edge.measure(); });
		return *std::max_element(lengths.begin(), lengths.end());
	}
};

using QuadrilateralElement2D = QuadrilateralElement<2>;
using QuadrilateralElement3D = QuadrilateralElement<3>;