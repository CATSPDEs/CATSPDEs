#pragma once
#include "Node.hpp"

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element

template <LocalIndex D, LocalIndex N>
using Element = std::array<Node<D>, N>; // abstract element

	// segments 

	template <LocalIndex D>
	using Segment = Element<D, 2>;

		using Segment1D = Segment<1>;
		using Segment2D = Segment<2>;
		using Segment3D = Segment<3>;

		inline double length(Segment1D const & s) {
			return fabs(s[1] - s[0]);
		}

		inline double length(Segment2D const & s) {
			return norm(s[1] - s[0]);
		}
	
	template <LocalIndex N>
	inline bool nodeInElement(Element<2, N> const & e, Node2D const & p) {
		// for convex 2D elements
		for (Index i = 0, j; i < N; ++i) {
			j = (i + 1) % N;
			if (crossProduct(e[i] - p, e[j] - p)[2] < 0.) return false;
		}
		return true;
	}

	// triangles

	template <LocalIndex D>
	using Triangle = Element<D, 3>;

		using Triangle2D = Triangle<2>;
		using Triangle3D = Triangle<3>;

		inline Node2D centroid(Triangle2D const & t) {
			return (t[0] + t[1] + t[2]) / 3.;
		}

		inline Triangle2D midNodes(Triangle2D const & t) {
			return { midNode(t[1], t[2]), midNode(t[0], t[2]), midNode(t[0], t[1]) };
		}

		inline std::array<Segment2D, 3> ribsOf(Triangle2D const & t) {
			return { { { t[1], t[2] }, { t[2], t[0] }, { t[0], t[1] } } };
		}

		inline double area(Triangle2D const & t) {
			return .5 * norm(crossProduct(t[1] - t[0], t[2] - t[0]));
		}

		inline double perimeter(Triangle2D const & t) {
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

	// quadrilaterals

	template <LocalIndex D>
	using Quadrilateral = Element<D, 4>;

		using Quadrilateral2D = Quadrilateral<2>;	
		using Quadrilateral3D = Quadrilateral<3>;