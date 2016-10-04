#include <unordered_map> // for red–green refinement
#include <algorithm>
#include <string>
#include "Triangulation.hpp"

Triangulation::Triangulation(Node const & lb, Node const & rt, double h) {
	/*
		author: 
			Alexander Žilyakov, May 2016
		edited:
			Oct 2016
		comments:
			here we construct dummy rect triangulation
			w/ left bottom point @lb and right top point @rt
			elements will be right-angled triangles
			w/ max hypotenuse size := @percent * min of rect sizes
			for @lb = origo, @rt = (3,2) and @percent = .9 
			it will look something like this:
					 ________(3,2)
					|\ |\ |\ |
					|_\|_\|_\|
					|\ |\ |\ |
					|_\|_\|_\|
				(0,0)  
	*/
	if (h <= 0) throw invalid_argument("max edge length must be a positive real number");
	Node size = rt - lb;
	// so x-projection of size is width of our rect and y-projection is height
	if (size.x() <= 0 || size.y() <= 0) throw invalid_argument("invalid rect");
	double hypotenuse = min(h, size.norm()); // hypotenuse and
	double dx = hypotenuse / sqrt(2.), // legs max sizes 
	       dy = dx;
	Index ix = ceil(size.x() / dx), // numb of INTERVALS on x-axis and
		   iy = ceil(size.y() / dy), // '' on y-axis, 
		   nx = ix + 1, // numb of NODES on x-axis and
		   ny = iy + 1, // on y-axis, 
		   i, j, t, // dummy indicies for looping over nodes and triangles
		   O, N, NW, E; // and here goes our so called STENCIL
	// O is center node, N is upper node (N stands for NORTH),
	// NW is north-west node, and E is east node of O
	// here we go:
	//
	//  NW __ N
	//    \   | \
	//     \  |  \
	//      \ |   \
	//        O __ E
	//
	bool canGoNorth, canGoNorther, canGoWest, canGoEast, canGoSouth; 
	// we will need this flags in order to check if we can add neighbor (of node or triangle)
	_nodes.resize(nx * ny); // so numb of points is nx * ny
	/* 
	_neighbors.resize(nx * ny - 1); // we will also construct vector of nodes neighbors here
	*/
	// just because it is straight-forward
	// in general we need to loop over triangles to construct it
	// recall that we will need neighbors in order to assemble portrait of CRS matrix
	// note that numb of elements here is nx * ny - 1 since 
	// if nodes i < j are neigbors, then _neighbors[i] contains j, but _neighbors[j] does not
	// because we have symmetry here! y'know, like, undirected graph
	// that is why we do not store neighbors of the very last node
	_triangles.resize(2 * ix * iy); // just draw and you see why
	dx = size.x() / ix; // normalizing (legs actual size)
	dy = size.y() / iy;
	hypotenuse = sqrt(dx * dx + dy * dy); // good old days w/ Pythagorean thm
	for (i = 0, O = 0, t = 0; i < ny; ++i)
		for (j = 0; j < nx; ++j, ++O) {
			_nodes[O] = lb + Node(j * dx, i * dy); // compute nodes
			canGoNorth = i + 1 < ny;
			canGoNorther = i + 2 < ny;
			canGoWest = j;
			canGoEast = j + 1 < nx;
			canGoSouth = i;
			N = O + nx; // north node!
			if (canGoNorth) { // add neighbors from the future!
				/*
				_neighbors[O].push_front(N);
				*/
				if (canGoWest) {
					/*
					_neighbors[O].push_front(NW = N - 1);
					*/
					NW = N - 1;
					_triangles[t].nodes(O, N, NW); // counterclockwise! you can check our STENCIL 
					// now we have to add neighbors of our triangle
					if (canGoNorther) _triangles[t].neighbors(0) = t + 2 * ix - 1; // just draw and you will see
					_triangles[t].neighbors(1) = t - 1;
					if (canGoEast) _triangles[t].neighbors(2) = t + 1;
					++t;
				}
			}
			if (canGoEast) {
				/*
				_neighbors[O].push_front(E = O + 1);
				*/
				E = O + 1;
				if (canGoNorth) {
					_triangles[t].nodes(O, E, N);
					_triangles[t].neighbors(0) = t + 1;
					if (canGoWest) _triangles[t].neighbors(1) = t - 1;
					if (canGoSouth) _triangles[t].neighbors(2) = t - 2 * ix + 1; // same--you have to draw it and look precisely on STENCIL
					++t;
				}
			}
		}
}

Triangulation::Triangulation(double h) 
	: _nodes({ Node(-1., 0.), Node(0., 0.), Node(0., 1.), Node(1., 0.), Node(0., -1.) })
	, _curves({ circleCurve })
	, _curvilinearEdges({
		CurvilinearEdge(0., .25, 0),
		CurvilinearEdge(.25, .5, 0),
		CurvilinearEdge(.5, .75, 0),
		CurvilinearEdge(.75, 1., 0)
	})
	, _triangles({
		Triangle(0, 1, 2, 3, -3, 1),
		Triangle(0, 4, 1, 2, 0, -4),
		Triangle(4, 3, 1, 3, 1, -5),
		Triangle(3, 2, 1, 0, 2, -2)
	}) {
	if (h < 1. / sqrt(2)) {
		if (h <= 0) throw invalid_argument("max edge length must be a positive real number");
		while (longestEdge() > h) refine();
	}
}

Triangulation::Triangulation(istream& nodes, istream& triangles, istream& neighbors) {
	Index n;
	nodes >> n; // numb of nodes
	_nodes.resize(n);
	nodes >> _nodes; // load nodes
	triangles >> n; // numb of triangles
	_triangles.resize(n);
	triangles >> _triangles;
	// fix neighbors
	for (Index i = 0; i < _triangles.size(); ++i)
		neighbors >> _triangles[i].neighbors();
}

Triangulation::Triangulation(istream& nodes, istream& triangles) {
	Index n;
	nodes >> n; // numb of nodes
	_nodes.resize(n);
	nodes >> _nodes; // load nodes
	triangles >> n; // numb of triangles
	_triangles.resize(n);
	triangles >> _triangles;
	// fix neighbors
	for (Index i = 0; i < _triangles.size(); ++i)
		for (Index j = i + 1; j < _triangles.size(); ++j)
			_makeNeighbors(i, j);
}

// private methods

bool Triangulation::_makeNeighbors(Index t1, Index t2) {
	// make _triangles[@t1] and ‘‘ @t2 neighbors
	// if they are not adjacent, return false
	array<array<LocalIndex, 2>, 2> commonNodes;
	LocalIndex i, j, k = 0;
	for (i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			if (_triangles[t1].nodes(i) == _triangles[t2].nodes(j)) {
				if (k > 1) // triangles share… more than 2 nodes?
					throw logic_error("invalid mesh: check out tringles #" + to_string(t1) + " and #" + to_string(t2));
				commonNodes[0][k] = i;
				commonNodes[1][k++] = j;
			}
	if (k < 2) return false; // triangles are not adjacent
	i = excludeIndicies(commonNodes[0][0], commonNodes[0][1]); // i := node of t1 against common edge 
	j = excludeIndicies(commonNodes[1][0], commonNodes[1][1]);
	_triangles[t1].neighbors(i) = t2;
	_triangles[t2].neighbors(j) = t1;
	return true; // we are gold
}

bool Triangulation::_checkNeighbor(Index t1, LocalIndex i, LocalIndex* j) {
	// check if ith neighbor of a t1th triangle also has _triangles[t1] as a neighbor
	// if so, j = local index of node of ith neighbor that is against ith node of t1th triangle
	SignedIndex t2 = _triangles[t1].neighbors(i);
	if (t2 < 0) return true;
	for (i = 0; i < 3; ++i)
		if (_triangles[t2].neighbors(i) == t1) {
			if (j) *j = i;
			return true;
		}
	return false;
}

SignedIndex Triangulation::_neighbor2edge(SignedIndex i) {
	// convert @i := _triangles[*].neighbor(**) to
	// index in _curvilinearEdges vector
	return - i - 2;
}

// public methods

AdjacencyList Triangulation::generateAdjList() {
	AdjacencyList adjList(numbOfNodes());
	array<Index, 3> elementNodesIndicies;
	for (Index i = 0; i < numbOfTriangles(); ++i) {
		elementNodesIndicies = l2g(i);
		sort(elementNodesIndicies.begin(), elementNodesIndicies.end());
		adjList[elementNodesIndicies[2]].insert(elementNodesIndicies[1]);
		adjList[elementNodesIndicies[2]].insert(elementNodesIndicies[0]);
		adjList[elementNodesIndicies[1]].insert(elementNodesIndicies[0]);
	}
	return adjList;
}

double Triangulation::area(Index i) { // compute area of ith triangle
	// we have norm of vector product underhood
	// neat!!
	Node u = _nodes[_triangles[i].nodes(1)] - _nodes[_triangles[i].nodes(0)];
	Node v = _nodes[_triangles[i].nodes(2)] - _nodes[_triangles[i].nodes(0)];
	return u.crossProductNorm(v) / 2;
}

Triangulation& Triangulation::save(ostream& outNodes, ostream& outTriangles) {
	outNodes.precision(15); // double precision
	outNodes << scientific << showpos;
	for (Node const & p : _nodes) outNodes << p;
	for (Triangle const & t : _triangles) outTriangles << t;
	return *this;
}

Triangulation& Triangulation::refine(Indicies& redList) {
	/*
		author: 
			Alexander Žilyakov, Jun 2016
		edited:
			Oct 2016
		comments:
			here we implement red–green refinement of our mesh
			@redList is a vector of indicies in _triangles to be red–refined
	*/
	unordered_map<Index, unsigned short> greenMap;
	// hash table of indicies of triangles w/ hanging nodes
	// we must refine them green later
	// key is index in _triangles vector,
	// value is numb of hanging nodes in (1, 2, or 3)
	Indicies redNeighborsList;
	// n <= 6 unknown neighbors of three new red triangles
	// that were added on previous iterations
	array<Index, 3> p, rp;
	// indicies of nodes of old triangle (i.e. triangle to be refined),
	// ‘‘ of new (red) nodes to be added,
	Index gp; // index of green node and green neighbor
	array<SignedIndex, 3> t, rt;
	// indicies of neighbor triangles of the old triangle,
	// ‘‘ of new (red) neighbor triangles
	Index i, j, k, m; // dummy indicies
	auto addExistingRedNodeFrom = [&](Index t) { // …from triangle _triangles[t]
		// we will need this function in order to 
		// organize red refinement if @redList contains neighbor triangles
		// because we do not want to add same nodes to _nodes vector several times
		for (m = 0; m < 3; ++m) if (_triangles[_triangles[t].neighbors(m)].neighbors(1) != i &&
									_triangles[_triangles[t].neighbors(m)].neighbors(2) != i) break;
		redNeighborsList.push_back(_triangles[t].neighbors(m + 1)); // add red neighbors from
		redNeighborsList.push_back(_triangles[t].neighbors(m + 2)); // previous refinements
		return _triangles[t].nodes(m);
	};
	// RED PART
	Indicies::iterator redListIter = redList.begin(),
					   searchStartIter = prev(redList.end()), 
					   searchIter;
	while (redListIter != redList.end()) {
		// well, since this iteration
		// _triangles[i] sholuld not ever be added to redList again
		greenMap[i = *redListIter] = 3; 
		// we will need ith nodes and neighbors later
		p = _triangles[i].nodes();
		t = _triangles[i].neighbors();
		// let’s deal w/ red points
		for (j = 0, k = 0; j < 3; ++j)
			if (!_checkNeighbor(i, j)) {
				// if our triangle has a neighbor which has already been red-refined,
				// then no need in creating a new node
				// we just have to find out its index
				rp[j] = addExistingRedNodeFrom(t[j]);
				++k; // how many nodes we do not need to create
			}
			else {
				// otherwise, well, let’s add a node!
				rp[j] = _nodes.size();
				if (t[j] < -1) { // curvilinear edges
					m = _neighbor2edge(t[j]);
					_nodes.push_back(_curves[_curvilinearEdges[m].curveIndex()](_curvilinearEdges[m].thetaMiddle()));
					_curvilinearEdges.push_back(CurvilinearEdge(_curvilinearEdges[m].thetaMiddle(), _curvilinearEdges[m].thetaEnd(), _curvilinearEdges[m].curveIndex()));
					_curvilinearEdges[m].thetaEnd() = _curvilinearEdges[m].thetaMiddle();
				}
				else _nodes.push_back(_nodes[p[excludeIndex(j)[0]]].midPoint(_nodes[p[excludeIndex(j)[1]]]));
			}
			if (k > 1 && (searchIter = find(next(searchStartIter), redList.end(), i)) != redList.end()) 
				redList.erase(searchIter); // if k equals 2 or 3
			// and we have also 3 red triangles yet to be added
			rt[2] = (rt[1] = (rt[0] = _triangles.size()) + 1) + 1;
			// our ith triangle splits into 4 new ones
			// central one will take place of the old one
			_triangles[i]
				.nodes(rp)
				// and its neighbors are known (3 other triangles) and will be added soon
				.neighbors(rt);
			// now we have to add 3 other triangles
			// you have to draw them not to get confused w/ numeration
			// or just TRUST ME I AM A DOCTOR
			_triangles.push_back(Triangle(p[0], rp[2], rp[1], i, t[1], t[2]));
			_triangles.push_back(Triangle(p[1], rp[0], rp[2], i, t[2], t[0]));
			_triangles.push_back(Triangle(p[2], rp[1], rp[0], i, t[0], t[1]));
			// fix curvilinear edges
			j = _curvilinearEdges.size();
			for (Index newTriangle : { _triangles.size() - 2, _triangles.size() - 3, _triangles.size() - 1 })
				if (_triangles[newTriangle].neighbors(1) < -1) 
					_triangles[newTriangle].neighbors(1) = _neighbor2edge(--j);
			// we set neighbors of our new triangles to be neighbors of 
			// our old (refined) triangle
			// so there n <= 6 neighbors to be found
			// if ith neighbors were red-refined previously 
			for (j = _triangles.size() - 3; j < _triangles.size(); ++j)
				for (Index neighbor : redNeighborsList)
					_makeNeighbors(j, neighbor);
			redNeighborsList.clear(); // we are done w/ neighbors
			// if ith has no neighbors or its neighbors were refined earlier, we are gold 
			// otherwise we have to deal w/ hanging nodes
			for (j = 0; j < 3; ++j)
				if (t[j] > -1 && ++greenMap[t[j]] == 2) 
						// if we have 2 (or 3 actually) hanging nodes,
						// we will not refine _triangle[t[j]] green
						redList.push_back(t[j]); // we will refine it red instead!
			// we do not want to search from the very begining
			// because O(n^2) is too slow
			// so we use this smart hack
			if (redListIter++ == searchStartIter) searchStartIter = prev(redList.end());
	}
	// GREEN PART
	for (auto const & keyValue : greenMap) {
		if (keyValue.second > 1) continue; // we should do green refinement iff there’s only one hanging node
		i = keyValue.first; // index of triangle to be green-refined
		for (j = 0; j < 3; ++j)
			if (!_checkNeighbor(i, j)) {
				rt[0] = _triangles[i].neighbors(j); // so _triangles[rt[0]] is red triangle w/ a handing node
				break;
			}
		gp = addExistingRedNodeFrom(rt[0]); // so _nodes[gp] is our hagning node
		// we want to split our ith triangle into two triangles
		// one of them we will add… 
		_triangles.push_back(Triangle(_triangles[i].nodes(j), _triangles[i].nodes(j + 1), gp,
									  -1, i, _triangles[i].neighbors(j + 2)));
		// …(fix green neighbor)…
		if (_triangles[i].neighbors(j + 2) > -1) 
			_makeNeighbors(_triangles.size() - 1, _triangles[i].neighbors(j + 2));
		// …and another one will take place of the old one
		_triangles[i]
			.nodes(_triangles[i].nodes(j), gp, _triangles[i].nodes(j + 2))
			.neighbors(-1, _triangles[i].neighbors(j + 1), _triangles.size() - 1);
		// finally, lets fix neighbors
		_makeNeighbors(i, redNeighborsList.front());
		_makeNeighbors(_triangles.size() - 1, redNeighborsList.back());
		redNeighborsList.clear();
	}
	return *this;
}

Triangulation& Triangulation::refine(unsigned numbOfRefinements) {
	/*
		author: 
			Alexander Žilyakov, Jun 2016
		edited:
			Oct 2016
		comments:
			 here we implement uniform refinement
			 works just like red green refinement, but 
			 ALL triangles will be refined in red fashion
			 so we do not need stuff like redList, greenMap etc. here
	*/
	Indicies redNeighborsList;
	array<Index, 3> p, rp;
	array<SignedIndex, 3> t, rt;
	Index i, j, m, n;
	auto addExistingRedNodeFrom = [&](Index t) {
		for (m = 0; m < 3; ++m) if (_triangles[_triangles[t].neighbors(m)].neighbors(1) != i &&
			_triangles[_triangles[t].neighbors(m)].neighbors(2) != i) break;
		redNeighborsList.push_back(_triangles[t].neighbors(m + 1));
		redNeighborsList.push_back(_triangles[t].neighbors(m + 2));
		return _triangles[t].nodes(m);
	};
	for (unsigned k = 0; k < numbOfRefinements; ++k) {
		n = _triangles.size();
		for (i = 0; i < n; ++i) {
			p = _triangles[i].nodes();
			t = _triangles[i].neighbors();
			for (j = 0; j < 3; ++j)
				if (!_checkNeighbor(i, j))
					rp[j] = addExistingRedNodeFrom(t[j]);
				else {
					rp[j] = _nodes.size();
					if (t[j] < -1) {
						m = _neighbor2edge(t[j]);
						_nodes.push_back(_curves[_curvilinearEdges[m].curveIndex()](_curvilinearEdges[m].thetaMiddle()));
						_curvilinearEdges.push_back(CurvilinearEdge(_curvilinearEdges[m].thetaMiddle(), _curvilinearEdges[m].thetaEnd(), _curvilinearEdges[m].curveIndex()));
						_curvilinearEdges[m].thetaEnd() = _curvilinearEdges[m].thetaMiddle();
					}
					else _nodes.push_back(_nodes[p[excludeIndex(j)[0]]].midPoint(_nodes[p[excludeIndex(j)[1]]]));
				}
				rt[2] = (rt[1] = (rt[0] = _triangles.size()) + 1) + 1;
				_triangles[i]
					.nodes(rp)
					.neighbors(rt);
				_triangles.push_back(Triangle(p[0], rp[2], rp[1], i, t[1], t[2]));
				_triangles.push_back(Triangle(p[1], rp[0], rp[2], i, t[2], t[0]));
				_triangles.push_back(Triangle(p[2], rp[1], rp[0], i, t[0], t[1]));
				j = _curvilinearEdges.size();
				for (Index newTriangle : { _triangles.size() - 2, _triangles.size() - 3, _triangles.size() - 1 })
					if (_triangles[newTriangle].neighbors(1) < -1)
						_triangles[newTriangle].neighbors(1) = _neighbor2edge(--j);
				for (j = _triangles.size() - 3; j < _triangles.size(); ++j)
					for (Index neighbor : redNeighborsList)
						_makeNeighbors(j, neighbor);
				redNeighborsList.clear();
		}
	}
	return *this;
}

vector<double> Triangulation::longestEdges() {
	vector<double> v(_triangles.size());
	for (Index i = 0; i < _triangles.size(); ++i) {
		v[i] = max(length(i, 0), length(i, 1));
		v[i] = max(v[i], length(i, 2));
	}
	return v;
}

double Triangulation::longestEdge() {
	vector<double> v = longestEdges();
	return *max_element(v.begin(), v.end());
}

vector<double> Triangulation::inscribedDiameters() {
	vector<double> v(_triangles.size());
	for (Index i = 0; i < _triangles.size(); ++i)
		v[i] = 4 * area(i) / perimeter(i);
	return v;
}

vector<double> Triangulation::qualityMeasure() {
	vector<double> v(_triangles.size()),
				   h(longestEdges()),
				   d(inscribedDiameters());
	for (Index i = 0; i < _triangles.size(); ++i)
		v[i] = sqrt(3) * d[i] / h[i]; 
	// we multiply by sqrt(3), so ideal triangles have quality measure = 1
	return v;
}

RibsNumeration Triangulation::computeRibsNumeration() {
	RibsNumeration res(_triangles.size());
	Index currentNumber = 0, n;
	LocalIndex k;
	for (Index i = 0; i < _triangles.size(); ++i)
		for (LocalIndex j = 0; j < 3; ++j) {
			n = _triangles[i].neighbors(j);
			// if neighbor of ith triangle has not been reached yet or does not exist, then we should numerate
			if (i < n || n < 0) res[i][j] = currentNumber++;
			// otherwise, rib already has a number so we need to use it
			else {
				if (!_checkNeighbor(i, j, &k)) throw logic_error("invalid mesh: check out tringles #" + to_string(i) + " and #" + to_string(n));
				res[i][j] = res[n][k];
			}
		}
	return res;
}

vector<Node> Triangulation::computeMiddleNodes(RibsNumeration const & ribs) {
	/*
		author:
			Alexander Žilyakov, Oct 2016
		edited:
			_	
		comments:
			we will not usually compute mid nodes in real FEM apps
			we will get them from elements if needed (e.g. in FEM apps
			using, for one, Nedelec or Crouzeix–Raviart finite elements,
			which DOFs numeration is induced by ribs’ / mid nodes’ numeration)

			anyway, it is useful to have this routine for the purpose of analysis
			we use it, for one, in Δ P1 CR — Δ P0 L Interpolant Visualization project
			in order to compute basis coefs u1Vec, u2Vec of velocity components u1, u2

			in real problems we will get such vectors solving a linear system
	*/
	vector<Node> res;
	for (Index i = 0; i < ribs.size(); ++i)
		for (LocalIndex j = 0; j < 3; ++j)
			if (ribs[i][j] + 1 > res.size()) res.push_back(
				_nodes[_triangles[i].nodes(j + 1)].midPoint(
					_nodes[_triangles[i].nodes(j + 2)]
				)
			);
	return res;
}

Boundary Triangulation::computeBoundary() {
	Boundary edges;
	for (Index i = 0; i < numbOfTriangles(); ++i) 
		for (LocalIndex edgeIndex : getBoundaryIndicies(i)) 
			edges.push_back({ _triangles[i].nodes(edgeIndex + 1), _triangles[i].nodes(edgeIndex + 2) });
	return edges;
}