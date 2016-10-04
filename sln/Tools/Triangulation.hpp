#pragma once
#include <iostream>
#include "Node.hpp"
#include "Triangle.hpp"
#include "Curve.hpp"
#include "CurvilinearEdge.hpp"
#include "AdjacencyList.hpp"
#include "vector.hpp"

/*
	Alexander Žilyakov, Aug 2016
*/

typedef vector<array<Index, 2>> Boundary;       
// we do not store boundary explicitly (since it is available from elements), but we can compute it if necessary
typedef vector<array<Index, 3>> RibsNumeration; 
// ditto for ribs; we might want to numerate ribs in order to implement: 
// * quadratic Lagrange triangles (i.e. we need numaration for midpoints of ribs),
// * any kind of FE w/ vector shape funcs (because they associated w/ ribs not nodes)
// RibsNumeration rn, rn[t][0] is the number of the rib against 0th node of tth triangle (see Triangulation def for details)
// such triangulation data structure (which stores ribs numeration but not ribs themselves) 
// is called “Nodes, Simple Ribs, and Triangles” 

class Triangulation { // this data structure is known as “Nodes and Triangles”
	// we have array of nodes := points on the plane and
	// array of triangles := array of nodes’ indicies + array of adjacent triangles’ indicies
	// we use because it is easy to:
	// * implement Delaunay algorithm [well, indeed]
	// * loop over elements (i.e. triangles) [stiffnes / mass matrix and load vector assembly]
	// * determine boundary edges [assembly of Robin BCs]
	// * refine mesh w/o reconstruction [we call it adaptive FEM and it is neat!] 
	// main elements:
	vector<Node> _nodes; // nodes (i.e. P-matrix),
	vector<Triangle> _triangles; // triangles (i.e. T-matrix or connectivity matrix)
	// for curvilinear refinement: 
	vector<Curve> _curves; // curves that make the boundary
	vector<CurvilinearEdge> _curvilinearEdges;
	// we will loop over our elements (i.e. over _triangles vector) to assemble stiffness matrix
	// no need to loop over boundary edges to assemble Robin BCs
	// because we can easily determine bndry while looping over elements (look at Triangle data structure!) 
	// in order to construct portrait of CRS(-like)-matrix, we also need to store neighbors of ith node 
	bool _checkNeighbor(Index, LocalIndex, LocalIndex* j = nullptr); // check if ith neighbor of a tth triangle also has _triangles[t] as a neighbor
	bool _makeNeighbors(Index, Index); // make 2 triangles neighbors
	SignedIndex _neighbor2edge(SignedIndex); // mapping between indicies
public:
	/* 
	(I) model domains constructors
	*/
	// (1) dummy rect triangulation
	Triangulation(Node const &, Node const &, double h = INFTY); 
	// (2) unit circle triangulation w/ center at origo,
	Triangulation(double h = INFTY); // @h := longest edge 
	// (3) import triangulation	
	// from nodes.dat, triangles.dat, and neighbors.dat
	Triangulation(istream&, istream&, istream&);
	// (4) ditto from nodes.dat and triangles.dat alone
	// O(m^2), m := numb of triangles (because we need to construct _neighbors list manually)
	Triangulation(istream&, istream&);
	/* 
	(II) simple inline methods
	*/
	double length(Index t, LocalIndex i) { // O(1)
		// compute length of ith edge of tth triangle
		return (_nodes[_triangles[t].nodes(i + 1)] - _nodes[_triangles[t].nodes(i + 2)]).norm();
	}
	Node getNode(Index i) const {
		return _nodes[i];
	}
	array<Node, 3> getNodes(Index t) { // ... of tth triangle
		return { _nodes[_triangles[t].nodes(0)], _nodes[_triangles[t].nodes(1)], _nodes[_triangles[t].nodes(2)] };
	}
	Node centroid(Index t) { // …of tth triangle
		return ( _nodes[_triangles[t].nodes(0)] + _nodes[_triangles[t].nodes(1)] + _nodes[_triangles[t].nodes(2)] ) / 3.;
	}
	array<Index, 3> l2g(Index t) const { // local to global nodes numeration
		return _triangles[t].nodes();
	}
	LocalIndicies getBoundaryIndicies(Index t) {
		LocalIndicies b;
		for (LocalIndex i = 0; i < 3; i++)
			if (_triangles[t].neighbors(i) < 0)
				b.push_back(i);
		return b;
	}
	double perimeter(Index t) { // O(1)
		return length(t, 0) + length(t, 1) + length(t, 2);
	}
	Index numbOfNodes() const { return _nodes.size(); }
	Index numbOfTriangles() const { return _triangles.size(); }
	/* 
	(III) more complex methods (exporting, refining etc.)
	*/
	AdjacencyList generateAdjList(); // nodes’ neighbors [we need it to construct portrait of CSlR matrix],
	double area(Index); // compute area of ith triangle
	Triangulation& save(ostream& nodes = cout, ostream& triangles = cout); // save mesh to std out
	Triangulation& refine(Indicies&); // red-green refinement
	Triangulation& refine(unsigned numbOfRefinements = 1); // uniform refinement
	/* 
	(IV) mesh quality measures
	*/
	// (1) compute vector of longest edges of all triangles
	vector<double> longestEdges(); // O(n), n := _triangles.size()
	double longestEdge();
	// (2) compute vector of diameters of inscribed circles of all triangles
	vector<double> inscribedDiameters(); // O(n), n := _triangles.size()
	vector<double> qualityMeasure(); // O(n)
	/*
	(VI) compute something that are not stored in our mesh explicitly
	*/
	RibsNumeration computeRibsNumeration(); // we may want to do this for some FE; see def of RibsNumeration 
	Boundary computeBoundary(); // O(numb of triangles); compute vector of arrays of indicies of _nodes that make the boudary 
	vector<Node> computeMiddleNodes(RibsNumeration const &);
	/* 
	(VII) user-defined mesh generator
	*/
	// useful for creating user-defined complex meshes
	// write your own generateMesh(), then use copy constructor:
	//     Triangulation Omega(generateMesh());
	friend Triangulation generateMesh();
};