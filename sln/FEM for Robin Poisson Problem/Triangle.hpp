#pragma once
#include <array>

// ssize_t
#if __BITS_PER_LONG != 64
typedef int ssize_t;
#else
typedef long ssize_t;
#endif

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
	size_t& nodes(size_t i) { return _nodes[i % 3]; }
	Triangle& nodes(size_t i, size_t j, size_t k) { 
		_nodes[0] = i;
		_nodes[1] = j;
		_nodes[2] = k;
		return *this; 
	}
	ssize_t& neighbors(size_t i) { return _neighbors[i % 3]; }
	Triangle& neighbors(ssize_t i, ssize_t j, ssize_t k) {
		_neighbors[0] = i;
		_neighbors[1] = j;
		_neighbors[2] = k;
		return *this;
	}
	ssize_t neighbor2node(size_t t) { // return index of the node against tth triangle
		for (size_t i = 0; i < 3; ++i)
			if (_neighbors[i] == t) return i;
		return -1; // well, t is not a neighbor
	}
	friend ostream& operator<<(ostream& out, Triangle const & t) {
		return out << t._nodes[0] << ' ' << t._nodes[1] << ' ' << t._nodes[2] << '\n';
	}
};