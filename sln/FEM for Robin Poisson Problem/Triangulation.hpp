#pragma once
#include <vector>
#include <forward_list>
#include <iostream>
#include "Node.hpp"
#include "Triangle.hpp"

typedef forward_list<size_t> Indicies;

class Triangulation { // this data structure is known as “Nodes and Triangles”
	// we have array of nodes := points on the plane and
	// array of triangles := array of nodes’ indicies + array of adjacent triangles’ indicies
	// we use because it is easy to:
	// * implement Delaunay algorithm [well, indeed]
	// * loop over elements (i.e. triangles) [stiffnes / mass matrix and load vector assembly]
	// * determine boundary edges [assembly of Robin BCs]
	// * refine mesh w/o reconstruction [we call it adaptive FEM and it is neat!] 
	vector<Node> _nodes; // nodes (i.e. P-matrix),
	vector<Indicies> _neighbors; // nodes’ neighbors [we need to construct it to assemble our CRS matrix],
	vector<Triangle> _triangles; // triangles (i.e. T-matrix or connectivity matrix), and
	double _h; // max size of triangle edge
	// we will loop over our elements (i.e. over _triangles vector) to assemble stiffness matrix
	// no need to loop over boundary edges to assemble Robin BCs
	// because we can easily determine bndry while looping over elements (look at Triangle data structure!) 
	// in order to construct portrait of CRS-matrix, we also need to store neighbors of ith node 
public:
	Triangulation(Node const &, Node const &, double percent = .5); // dummy rect triangulation
	double length(size_t t, size_t i) { // compute length of ith edge of tth triangle
		return (_nodes[_triangles[t].nodes((i + 1) % 3)] - _nodes[_triangles[t].nodes((i + 2) % 3)]).norm();
	}
	double area(size_t); // compute area of ith triangle
	Triangulation& save(ostream& nodes = cout, ostream& triangles = cout); // save mesh to std out
	Triangulation& refine(Indicies&); // red-green refinement
};