#pragma once
#include <array>

// ssize_t
#if __BITS_PER_LONG != 64
typedef int ssize_t;
#else
typedef long ssize_t;
#endif

typedef unsigned short localIndex; // in {0, 1, 2}

using namespace std;

class Triangle {
	array<size_t, 3> _nodes; // nodes in triangle, counterclockwise! 
	array<ssize_t, 3> _neighbors; // neighbors of the triangle (they are... triangles!)
	// ith triangle against ith node in _nodes array
	// if there is no triangle (i.e. edge is part of boundary), then index is -1
	// that is why we use ssize_t instead of size_t here
public:
	explicit Triangle(size_t p1 = 0, size_t p2 = 0, size_t p3 = 0,
					  ssize_t t1 = -1, ssize_t t2 = -1, ssize_t t3 = -1) 
		: _nodes{ {p1, p2, p3} }
		, _neighbors{ {t1, t2, t3} } {}
	array<size_t, 3> nodes() const { return _nodes; }
	Triangle& nodes(array<size_t, 3> const & newNodes) {
		_nodes = newNodes;
		return *this;
	}
	size_t& nodes(localIndex i) { return _nodes[i % 3]; }
	Triangle& nodes(size_t i, size_t j, size_t k) {
		_nodes[0] = i;
		_nodes[1] = j;
		_nodes[2] = k;
		return *this;
	}
	array<ssize_t, 3> neighbors() { return _neighbors; }
	Triangle& neighbors(array<ssize_t, 3> const & newNeighbors) {
		_neighbors = newNeighbors;
		return *this;
	}
	ssize_t& neighbors(localIndex i) { return _neighbors[i % 3]; }
	Triangle& neighbors(ssize_t i, ssize_t j, ssize_t k) {
		_neighbors[0] = i;
		_neighbors[1] = j;
		_neighbors[2] = k;
		return *this;
	}
	friend ostream& operator<<(ostream& out, Triangle const & t) {
		return out << t._nodes[0] << ' ' << t._nodes[1] << ' ' << t._nodes[2] << '\n';
	}
};

inline array<localIndex, 2> excludeIndex(localIndex i) { 
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

inline localIndex excludeIndicies(localIndex i, localIndex j) { 
	// exclude i, j from <0, 1, 2>		
	if (i > j) swap(i, j);
	if (i == 0) return j == 1 ? 2 : 1;
	return 0;
}