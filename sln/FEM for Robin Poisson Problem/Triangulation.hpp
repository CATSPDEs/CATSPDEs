#pragma once
#include <vector>
#include <list>
#include <iostream>
#include "Node.hpp"
#include "Triangle.hpp"
#include "Curve.hpp"
#include "CurvilinearEdge.hpp"

typedef list<size_t> Indicies;
typedef list<unsigned short> LocalIndicies;

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
	vector<Curve> _curves; // curves that make the boundary
	vector<CurvilinearEdge> _edges;
	// we will loop over our elements (i.e. over _triangles vector) to assemble stiffness matrix
	// no need to loop over boundary edges to assemble Robin BCs
	// because we can easily determine bndry while looping over elements (look at Triangle data structure!) 
	// in order to construct portrait of CRS-matrix, we also need to store neighbors of ith node 
	bool _checkNeighbor(size_t, localIndex); // check if ith neighbor of a tth triangle also has _triangles[t] as a neighbor
	bool _makeNeighbors(size_t, size_t); // make 2 triangles neighbors
	ssize_t _neighbor2edge(ssize_t); // mapping between indicies
	// mesh quality measures
	vector<double> _longestEdges(); // O(n), n := _triangles.size()
	// compute vector of longest edges of all triangles
	vector<double> _inscribedDiameters(); // O(n), n := _triangles.size()
	// compute vector of diameters of inscribed circles of all triangles
public:
	Triangulation(Node const &, Node const &, double percent = .5); // dummy rect triangulation
	Triangulation(vector<Node> const &, vector<Triangle> const &, 
				  vector<Curve> const &, vector<CurvilinearEdge> const &);
	double length(size_t t, localIndex i) { // O(1)
		// compute length of ith edge of tth triangle
		return (_nodes[_triangles[t].nodes(i + 1)] - _nodes[_triangles[t].nodes(i + 2)]).norm();
	}
	array<Node, 3> getNodes(size_t t) {
		return { _nodes[_triangles[t].nodes(0)],_nodes[_triangles[t].nodes(1)],_nodes[_triangles[t].nodes(2)] };
	}
	array<size_t, 3> l2g(size_t t) { // local to global nodes numeration
		return _triangles[t].nodes();
	}
	LocalIndicies getBoundaryIndicies(size_t t) {
		LocalIndicies b;
		for (localIndex i = 0; i < 3; i++)
			if (_triangles[t].neighbors(i) < 0)
				b.push_back(i);
		return b;
	}
	double perimeter(size_t t) { // O(1)
		return length(t, 0) + length(t, 1) + length(t, 2);
	}
	size_t numbOfNodes() const { return _nodes.size(); }
	size_t numbOfTriangles() const { return _triangles.size(); }
	double area(size_t); // compute area of ith triangle
	Triangulation& save(ostream& nodes = cout, ostream& triangles = cout); // save mesh to std out
	Triangulation& refine(Indicies&); // red-green refinement
	Triangulation& refine(unsigned numbOfRefinements = 1); // uniform refinement
	vector<double> qualityMeasure(); // O(n)
};

inline localIndex nextIndex(localIndex i) {
	return (i + 1) % 3;
}