#pragma once
#include "Node.hpp"

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element

template <LocalIndex D, LocalIndex N>
using Element = std::array<Node<D>, N>; // abstract element

using Triangle = Element<2, 3>;

	inline Node2D centroid(Triangle const & t) {
		return (t[0] + t[1] + t[2]) / 3.;
	}

	inline Triangle midNodes(Triangle const & t) {
		return { midNode(t[1], t[2]), midNode(t[0], t[2]), midNode(t[0], t[1]) };
	}

	inline double area(Triangle const & t) {
		return .5 * norm(crossProduct(t[1] - t[0], t[2] - t[0]));
	}

	inline double perimeter(Triangle const & t) {
		return norm(t[1] - t[0]) + norm(t[2] - t[1]) + norm(t[0] - t[2]);
	}

	inline std::array<LocalIndex, 2> excludeIndex(LocalIndex i) {
		// exclude i from <0, 1, 2> counterclockwise		
		//   _____<_____
		//  |     2	    |
		//  |   /   \	|
		//  |  0 ___ 1  |
		//  |_____>_____|
		//
		i %= 3;
		if (i == 0) return { 1, 2 };
		if (i == 1) return { 2, 0 };
		return { 0, 1 };
	}

	inline LocalIndex excludeIndicies(LocalIndex i, LocalIndex j) {
		// exclude i, j from <0, 1, 2>		
		if (i > j) std::swap(i, j);
		if (i == 0) return j == 1 ? 2 : 1;
		return 0;
	}

	inline LocalIndex nextIndex(LocalIndex i) {
		return (i + 1) % 3;
	}
