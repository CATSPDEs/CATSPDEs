#pragma once
#include <numeric> // inner_product
#include "array.hpp"

template <LocalIndex N> 
using Node = std::array<double, N>; // abstract node

	// dot product
	template <LocalIndex N>
	double operator*(Node<N> const & u, Node<N> const & v) {
		return std::inner_product(u.begin(), u.end(), v.begin(), 0.);
	}

	// norm
	template <LocalIndex N>
	double norm(Node<N> const & u) {
		return sqrt(u * u);
	}

	// middle point
	template <LocalIndex N>
	Node<N> midNode(Node<N> const & u, Node<N> const & v) {
		return .5 * (u + v);
	}

using Node1D = Node<1>;

using Node2D = Node<2>;

using Node3D = Node<3>;

	inline Node3D crossProduct(Node3D const & u, Node3D const & v) {
		return {
			u[1] * v[2] - u[2] * v[1],
			u[2] * v[0] - u[0] * v[2],
			u[0] * v[1] - u[1] * v[0]
		};
	}

	inline Node3D crossProduct(Node2D const & u, Node2D const & v) {
		return { 0., 0, u[0] * v[1] - u[1] * v[0] };
	}