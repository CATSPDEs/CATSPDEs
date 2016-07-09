#pragma once
#include "array.hpp"
#include "indicies.hpp"

class Triangle {
	array<Index, 3> _nodes; // nodes in triangle, counterclockwise! 
	array<SignedIndex, 3> _neighbors; // neighbors of the triangle (they are... triangles!)
	// ith triangle against ith node in _nodes array
	// if there is no triangle (i.e. edge is part of boundary), then index is < 0
	// that is why we use SignedIndex instead of Index here
public:
	explicit Triangle(Index p1 = 0, Index p2 = 0, Index p3 = 0,
					  SignedIndex t1 = -1, SignedIndex t2 = -1, SignedIndex t3 = -1) 
		: _nodes{ {p1, p2, p3} }
		, _neighbors{ {t1, t2, t3} } {}
	// (1) nodes methods
	array<Index, 3>& nodes() { return _nodes; } // set / get
	array<Index, 3>  nodes() const { return _nodes; } // get
	Triangle&        nodes(array<Index, 3> const & newNodes) { // set
		_nodes = newNodes;
		return *this;
	}
	Index&           nodes(LocalIndex i) { return _nodes[i % 3]; } // set / get
	Triangle&        nodes(Index i, Index j, Index k) { // set
		_nodes[0] = i;
		_nodes[1] = j;
		_nodes[2] = k;
		return *this;
	}
	// (2) neighbors methods
	array<SignedIndex, 3>& neighbors() { return _neighbors; } // set / get
	array<SignedIndex, 3>  neighbors() const { return _neighbors; } // get
	Triangle&              neighbors(array<SignedIndex, 3> const & newNeighbors) { // set
		_neighbors = newNeighbors;
		return *this;
	}
	SignedIndex&           neighbors(LocalIndex i) { return _neighbors[i % 3]; } // set / get
	Triangle&              neighbors(SignedIndex i, SignedIndex j, SignedIndex k) { // set
		_neighbors[0] = i;
		_neighbors[1] = j;
		_neighbors[2] = k;
		return *this;
	}
	// (3) import / export
	friend ostream& operator<<(ostream& output, Triangle const & t) {
		return output << t._nodes[0] << ' ' << t._nodes[1] << ' ' << t._nodes[2] << '\n';
	}
	friend istream& operator>>(istream& input, Triangle& t) {
		return input >> t._nodes[0] >> t._nodes[1] >> t._nodes[2];
	}
};

inline array<LocalIndex, 2> excludeIndex(LocalIndex i) { 
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
	else return { 0, 1 };
}

inline LocalIndex excludeIndicies(LocalIndex i, LocalIndex j) { 
	// exclude i, j from <0, 1, 2>		
	if (i > j) swap(i, j);
	if (i == 0) return j == 1 ? 2 : 1;
	return 0;
}

inline LocalIndex nextIndex(LocalIndex i) {
	return (i + 1) % 3;
}