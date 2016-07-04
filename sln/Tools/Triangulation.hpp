#pragma once
#include <list>
#include <iostream>
#include "Node.hpp"
#include "Triangle.hpp"
#include "Curve.hpp"
#include "CurvilinearEdge.hpp"
#include "AdjacencyList.hpp"

double const INFTY = 10e50;

typedef list<size_t> Indicies;
typedef list<unsigned short> LocalIndicies;
typedef vector<array<size_t, 2>> Boundary; // we do not save boudary explicitly (since it is available from elements), but we can compute it if necessary

class Triangulation { // this data structure is known as “Nodes and Triangles”
	// we have array of nodes := points on the plane and
	// array of triangles := array of nodes’ indicies + array of adjacent triangles’ indicies
	// we use because it is easy to:
	// * implement Delaunay algorithm [well, indeed]
	// * loop over elements (i.e. triangles) [stiffnes / mass matrix and load vector assembly]
	// * determine boundary edges [assembly of Robin BCs]
	// * refine mesh w/o reconstruction [we call it adaptive FEM and it is neat!] 
	vector<Node> _nodes; // nodes (i.e. P-matrix),
	vector<Triangle> _triangles; // triangles (i.e. T-matrix or connectivity matrix), and
	vector<Curve> _curves; // curves that make the boundary
	vector<CurvilinearEdge> _edges;
	// we will loop over our elements (i.e. over _triangles vector) to assemble stiffness matrix
	// no need to loop over boundary edges to assemble Robin BCs
	// because we can easily determine bndry while looping over elements (look at Triangle data structure!) 
	// in order to construct portrait of CRS-matrix, we also need to store neighbors of ith node 
	bool _checkNeighbor(size_t, localIndex); // check if ith neighbor of a tth triangle also has _triangles[t] as a neighbor
	bool _makeNeighbors(size_t, size_t); // make 2 triangles neighbors
	ssize_t _neighbor2edge(ssize_t); // mapping between indicies
public:
	/* (I) model domains constructors
	*/
	// (1) dummy rect triangulation
	Triangulation(Node const &, Node const &, double percent = .5); 
	// (2) unit circle triangulation w/ center at origo,
	Triangulation(double h = INFTY); // @h := longest edge 
	Triangulation(vector<Node> const &, vector<Triangle> const &, 
				  vector<Curve> const &, vector<CurvilinearEdge> const &); // manually constructed mesh
	Triangulation(istream&, istream&); // import triangulation
	// from nodes.dat and triangles.dat
	// O(m^2), m := numb of triangles (because we need to construct _neighbors list)
	double length(size_t t, localIndex i) { // O(1)
		// compute length of ith edge of tth triangle
		return (_nodes[_triangles[t].nodes(i + 1)] - _nodes[_triangles[t].nodes(i + 2)]).norm();
	}
	/* (II) simple inline methods
	*/
	Node getNode(size_t i) const {
		return _nodes[i];
	}
	array<Node, 3> getNodes(size_t t) { // ... of tth triangle
		return { _nodes[_triangles[t].nodes(0)],_nodes[_triangles[t].nodes(1)],_nodes[_triangles[t].nodes(2)] };
	}
	array<size_t, 3> l2g(size_t t) const { // local to global nodes numeration
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
	/* (III) more complex methods (exporting, refining etc.)
	*/
	AdjacencyList generateAdjList(); // nodes’ neighbors [we need it to construct portrait of CSlR matrix],
	double area(size_t); // compute area of ith triangle
	Triangulation& save(ostream& nodes = cout, ostream& triangles = cout); // save mesh to std out
	Triangulation& refine(Indicies&); // red-green refinement
	Triangulation& refine(unsigned numbOfRefinements = 1); // uniform refinement
	/* (IV) mesh quality measures
	*/
	// (1) compute vector of longest edges of all triangles
	vector<double> longestEdges(); // O(n), n := _triangles.size()
	double longestEdge();
	// (2) compute vector of diameters of inscribed circles of all triangles
	vector<double> inscribedDiameters(); // O(n), n := _triangles.size()
	vector<double> qualityMeasure(); // O(n)
	Boundary computeBoundary(); // O(numb of triangles); compute vector of arrays of indicies of _nodes that make the boudary 
	/* (VI) user-defined mesh generator
	*/
	// useful for creating user-defined complex meshes
	// write your own generateMesh(), then use copy constructor:
	//     Triangulation Omega(generateMesh());
	friend Triangulation generateMesh();
};

inline localIndex nextIndex(localIndex i) {
	return (i + 1) % 3;
}