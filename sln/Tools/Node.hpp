#pragma once
#include <numeric> // inner_product
#include "array.hpp"

// D := dimension of the node (1 for segments, 2 for triangulations etc.)
// abstract node
template <LocalIndex D> 
using AbstractNode = std::array<double, D>;

	// dot product
	template <LocalIndex D>
	double operator*(AbstractNode<D> const & u, AbstractNode<D> const & v) {
		return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
	}

	// norm
	template <LocalIndex D>
	double norm(AbstractNode<D> const & u) {
		return sqrt(u * u);
	}

	// middle point
	template <LocalIndex D>
	AbstractNode<D> midNode(AbstractNode<D> const & u, AbstractNode<D> const & v) {
		return .5 * (u + v);
	}

	// ×–product (3D)
	inline AbstractNode<3> crossProduct(AbstractNode<3> const & u, AbstractNode<3> const & v) {
		return {
			u[1] * v[2] - u[2] * v[1],
			u[2] * v[0] - u[0] * v[2],
			u[0] * v[1] - u[1] * v[0]
		};
	}

	// ×–product (2D)
	inline AbstractNode<3> crossProduct(AbstractNode<2> const & u, AbstractNode<2> const & v) {
		return { 0., 0, u[0] * v[1] - u[1] * v[0] };
	}

// we want Node<1> to be just double, not array<double, 1>
// so here goes some magic
template <LocalIndex D> struct NodeTypedef { typedef AbstractNode<D> type; };
template<> struct NodeTypedef<1> { typedef double type; };
// abstract node
template <LocalIndex D>
using Node = typename NodeTypedef<D>::type;

using Node2D = Node<2>;
using Node3D = Node<3>;