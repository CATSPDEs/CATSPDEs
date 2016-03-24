#pragma once
#include <array>
#include <forward_list>
#include <vector>
#include "Node.hpp"

using namespace std;
typedef array<size_t, 3> Triangle;
typedef array<size_t, 2> Edge;
typedef forward_list<size_t> Neighbors;

class Triandulation { // our mesh is 
	vector<Node> _nodes; // nodes (i.e. P-matrix),
	vector<Neighbors> _neighbors; // nodes’ neighbors,
	vector<Triangle> _triangles; // triangles (i.e. T-matrix or connectivity matrix), and
	vector<Edge> _boundary; // boundary edges
	double _h; // max size of triangle edge
	// we will loop over our elements (i.e. over _triangles vector) to assembly stiffness matrix
	// and over boundary edges to assembly Robin BCs
	// so it is useful to store boundary edges
	// in order to construct portrait of CRS-matrix, we also need to store neighbors of ith node 
public:
	Triandulation(Node const & lb, Node const & rt, double percent = .5) {
		if (1 <= percent || percent <= 0) throw invalid_argument("3rd parameter should be el of (0, 1)");
		Node size = rt - lb;
		if (!size.x() || !size.y()) throw invalid_argument("invalid rect");
		_h = percent * min(size.x(), size.y()); // hypotenuse and
		double dx = _h / sqrt(2), // legs max sizes 
			   dy = dx;
		size_t nx = ceil(size.x() / dx) + 1, // numb of nodes on x-axis and
			   ny = ceil(size.y() / dy) + 1, i, j, O, N, NW, E; // on y-axis, 
		_nodes.resize(nx * ny); // so numb of points = nx * ny
		_neighbors.resize(nx * ny);
		_triangles.resize(2 * (nx - 1) * (ny - 1));
		dx = size.x() / (nx - 1); // normalizing (legs actual size)
		dy = size.y() / (ny - 1);
		// auto north = [&](size_t i) { return i + nx; };
		// todo: stencils!!
		for (i = 0, O = 0; i < ny; ++i)
			for (j = 0; j < nx; ++j, ++O) {
				_nodes[O] = lb + Node(j * dx, i * dy); // compute nodes
				N = O + nx;
				if (i + 1 < ny) { // add neighbors from the future!
					_neighbors[O].push_front(N);
					if (j) {
						_neighbors[O].push_front(NW = N - 1);
						_triangles.push_back({O, N, NW}); // counterclockwise
					}
				}
				if (j + 1 < nx) {
					_neighbors[O].push_front(E = O + 1);
					if (i + 1 < ny) _triangles.push_back({ O, E, N });
				}
			}
	}
	double length(size_t i) { // compute length of ith edge
		return (_nodes[_boundary[i][1]] - _nodes[_boundary[i][0]]).norm();
	}
	double area(size_t i) { // compute area of ith triangle; norm of vector product underhood
		Node u = _nodes[_triangles[i][1]] - _nodes[_triangles[i][0]];
		Node v = _nodes[_triangles[i][2]] - _nodes[_triangles[i][0]];
		return u.crossProductNorm(v) / 2;
	}
};