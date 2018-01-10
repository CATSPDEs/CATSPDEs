#pragma once
#include "boost/optional.hpp"
#include "Mapping.hpp"
#include "AbstractMesh.hpp"
#include "CurvilinearEdge.hpp"
#include "SymmetricCSlCMatrix.hpp" // for ribs numeration

/*
	Alexander Žilyakov, Aug 2016
*/

class Triangulation 
	: public AbstractMesh<2, 3> { 
	// this data structure is known as “Nodes and Triangles”
	// we have array of nodes := points on the plane and
	// array of triangles := array of nodes’ indicies + array of adjacent triangles’ indicies
	// we use because it is easy to:
	// * implement Delaunay algorithm [well, indeed]
	// * loop over elements (i.e. triangles) [stiffnes / mass matrix and load vector assembly]
	// * determine boundary edges [assembly of Robin BCs]
	// * refine mesh w/o reconstruction [we call it adaptive FEM and it is neat!] 	
	std::vector<std::array<SignedIndex, 3>> _neighbors;

	std::pair<Index, std::vector<std::array<Index, 3>>> _ribs;

	SymmetricCSlCMatrix<double>* _ribsNumeration = nullptr;
	
	// for curvilinear refinement: 
	std::vector<Curve2D> _curves; // curves that make the boundary
	std::vector<CurvilinearEdge> _edges;
	// we will loop over our elements (i.e. over _triangles vector) to assemble stiffness matrix
	// no need to loop over boundary edges to assemble Robin BCs
	// because we can easily determine bndry while looping over elements (look at Triangle data structure!) 
	// in order to construct portrait of CRS(-like)-matrix, we also need to store neighbors of ith node 
	std::vector<std::vector<Index>> _fineNeighborsIndicies;
	bool _makeNeighbors(Index, Index); // make 2 triangles neighbors
	SignedIndex _neighbor2edge(SignedIndex) const; // mapping between indicies
public:
	~Triangulation() { if (_ribsNumeration) delete _ribsNumeration; }
	Triangulation() {}
	Triangulation(
		std::vector<Node2D> const & nodes,
		std::vector<std::array<Index, 3>> const & elements,
		std::vector<std::array<SignedIndex, 3>> const & neighbors,
		std::vector<Curve2D> const & curves,
		std::vector<CurvilinearEdge> const & edges
	) 
		: AbstractMesh(nodes, elements)
		, _fineNeighborsIndicies(elements.size())
		, _neighbors(neighbors)
		, _curves(curves)
		, _edges(edges)
	{}
	Triangulation(Node2D const &, Node2D const &, Index ndx = 1, Index ndy = 1);
	// get numb of ribs 
	Index numbOfRibs() const {
		return _ribsNumeration->nnz() - numbOfNodes();
		// return _ribs.first;
	}
	// get indicies of ribs of tth triangle
	std::array<Index, 3> getRibsIndicies(Index t) const {
		auto nodesIndicies = getNodesIndicies(t);
		return { 
			_ribsNumeration->lvalIndex(nodesIndicies[1], nodesIndicies[2]),
			_ribsNumeration->lvalIndex(nodesIndicies[0], nodesIndicies[2]),
			_ribsNumeration->lvalIndex(nodesIndicies[0], nodesIndicies[1])
		};
		//return _ribs.second[t];
	}
	// get ribs numeration
	std::vector<std::array<Index, 3>> getRibsNumeration() const {
		std::vector<std::array<Index, 3>> res(numbOfRibs());
		for (Index t = 0; t < numbOfElements(); ++t) res[t] = getRibsIndicies(t);
		return res;
		//return _ribs.second;
	}
	// get indicies of neighbors of tth triangle
	std::array<SignedIndex, 3> getNeighborsIndicies(Index t) const {
		return _neighbors[t];
	}
	// get local rib indicies (if triangles are adjacent)
	boost::optional<std::array<LocalIndex, 2>> Triangulation::getCommonRibLocalIndicies(Index t1, Index t2) const;
	// get vector of local indicies of nodes triangles t1 and t2 share
	std::vector<std::array<LocalIndex, 2>> Triangulation::getCommonNodesLocalIndicies(Index t1, Index t2) const;
	// virtual methods to be implemented
	std::vector<Index> getFineNeighborsIndicies(Index e) const final {
		return _fineNeighborsIndicies[e];
	}
	Triangulation& import(std::istream& from = std::cin) final;
	void export(std::ostream& to = std::cout, Parameters const & params = {}) const final;
	// in order to work w/ strings, not streams
	using AbstractMesh::import;
	using AbstractMesh::export;
	Triangulation& refine(Index numbOfRefinements = 1) final;

	Triangulation& refine(Indicies&); // red—green refinement
	Triangulation& computeNeighbors(); // O(m^2), m := numb of triangles (because we need to construct _neighbors list manually)
	Triangulation& enumerateRibs(); // O(m), m := numb of triangles

	Node2D getRibNode(Index, LocalIndex, double) const;
	
	Triangulation& truncate(Predicate2D const p) {
		for (Index e = 0; e < numbOfElements(); ++e)
			if (!p(centroid(getElement(e)))) _ghostElements.insert(e);
		if (!_neighbors.size()) computeNeighbors();
		for (Index e = 0; e < numbOfElements(); ++e)
			for (auto& n : _neighbors[e])
				if (n >= 0 && _ghostElements.find(n) != _ghostElements.end()) n = -1;
		return *this;
	}

	/* 
	(I) model domains constructors
	*/
	// (1) dummy rect triangulation
	//Triangulation(Node const &, Node const &, double h = INFTY); 
	// (2) unit circle triangulation w/ center at origo,
	//Triangulation(double h = INFTY); // @h := longest edge 
	// (3) import triangulation	
	// from nodes.dat, triangles.dat, and neighbors.dat
	//Triangulation(istream&, istream&, istream&);
	// (4) ditto from nodes.dat and triangles.dat alone
	// O(m^2), m := numb of triangles (because we need to construct _neighbors list manually)
	//Triangulation(istream&, istream&);

	/* 
	(II) simple inline methods
	*/
	//double perimeter(Index t) { // O(1)
	//	return length(t, 0) + length(t, 1) + length(t, 2);
	//}

	
	/* 
	(III) more complex methods (exporting, refining etc.)
	*/
	//AdjacencyList generateAdjList(); // nodes’ neighbors [we need it to construct portrait of CSlR matrix],
	//Triangulation& save(ostream& nodes = cout, ostream& triangles = cout); // save mesh to std out
	//Triangulation& refine(Indicies&); // red-green refinement
	//Triangulation& refine(unsigned numbOfRefinements = 1); // uniform refinement
	//array<LocalIndex, 2> Triangulation::getCommonRibLocalIndicies (Index, Index) const;
	//array<LocalIndex, 2> Triangulation::getCommonNodeLocalIndicies(Index, Index) const;

	/* 
	(IV) mesh quality measures
	*/
	// (1) compute vector of longest edges of all triangles
	//vector<double> longestEdges(); // O(n), n := _triangles.size()
	//double longestEdge();
	// (2) compute vector of diameters of inscribed circles of all triangles
	//vector<double> inscribedDiameters(); // O(n), n := _triangles.size()
	//vector<double> qualityMeasure(); // O(n)

	/*
	(VI) compute something that are not stored in our mesh explicitly
	*/
	//RibsNumeration computeRibsNumeration(Index* numbOfRibs = nullptr, Index* numbOfBndryRibs = nullptr); // we may want to do this for some FE; see def of RibsNumeration 
	//Boundary computeBoundary(); // O(numb of triangles); compute vector of arrays of indicies of _nodes that make the boudary 
	//vector<Node> computeMiddleNodes(RibsNumeration const &);
	
	/* 
	(VII) user-defined mesh generator
	*/
	// useful for creating user-defined complex meshes
	// write your own generateMesh(), then use copy constructor:
	//     Triangulation Omega(generateMesh());
	//friend Triangulation generateMesh();
};