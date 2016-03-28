#include <algorithm> // find()
#include "Triangulation.hpp"

Triangulation::Triangulation(Node const & lb, Node const & rt, double percent) {
	// here we construct dummy rect triangulation
	// w/ left bottom point @lb and right top point @rt
	// elements will be right-angled triangles
	// w/ max hypotenuse size := @percent * min of rect sizes
	// for @lb = origo, @rt = (3,2) and @percent = .9 
	// it will look something like this:
	//   ________(3,2)
	//  |\ |\ |\ |
	//	|_\|_\|_\|
	//  |\ |\ |\ |
	//	|_\|_\|_\|
	//  (0,0)    
	if (1 <= percent || percent <= 0) throw invalid_argument("3rd parameter should be el of (0, 1)");
	Node size = rt - lb;
	// so x-projection of size is width of our rect and y-projection is height
	if (size.x() <= 0 || size.y() <= 0) throw invalid_argument("invalid rect");
	_h = percent * min(size.x(), size.y()); // hypotenuse and
	double dx = _h / sqrt(2), // legs max sizes 
		dy = dx;
	size_t ix = ceil(size.x() / dx), // numb of INTERVALS on x-axis and
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
	_neighbors.resize(nx * ny - 1); // we will also construct vector of nodes neighbors here
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
	_h = sqrt(dx * dx + dy * dy); // good old days w/ Pythagorean thm
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
				_neighbors[O].push_front(N);
				if (canGoWest) {
					_neighbors[O].push_front(NW = N - 1);
					_triangles[t].nodes(O, N, NW); // counterclockwise! you can check our STENCIL 
					// now we have to add neighbors of our triangle
					if (canGoNorther) _triangles[t].neighbors(0) = t + 2 * ix - 1; // just draw and you will see
					_triangles[t].neighbors(1) = t - 1;
					if (canGoEast) _triangles[t].neighbors(2) = t + 1;
					++t;
				}
			}
			if (canGoEast) {
				_neighbors[O].push_front(E = O + 1);
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

double Triangulation::area(size_t i) { // compute area of ith triangle
	// we have norm of vector product underhood
	// neat!!
	Node u = _nodes[_triangles[i].nodes(1)] - _nodes[_triangles[i].nodes(0)];
	Node v = _nodes[_triangles[i].nodes(2)] - _nodes[_triangles[i].nodes(0)];
	return u.crossProductNorm(v) / 2;
}

Triangulation& Triangulation::save(ostream& outNodes, ostream& outTriangles) {
	outNodes.precision(15); // double precision
	outNodes << scientific;
	for (Node p : _nodes) outNodes << p;
	for (Triangle t : _triangles) outTriangles << t;
	return *this;
}

Triangulation& Triangulation::refine(Indicies& indicies) {
	// @indicies is vector of indicies in _triangles to be red-refined

	Indicies::iterator iter1, iter2;
	size_t i, 
		   p1, p2, p3, 
		   t1, t2, t3,
		   rp1, rp2, rp3,
		   rt1, rt2, rt3,
		   gpl1, gpl2, gpl3,
		   gtn; // dummy indicies
	Triangle old;

	auto green = [&](size_t t, size_t first, size_t splittingNode, size_t second) { // green refinement of tth triangle
		// w/ @first and @second new neighbors
		// and new node @splittingNode
		gpl3 = (gpl2 = (gpl1 = _triangles[t].neighbor2node(i)) + 1) + 1;
		// gpl1 := first LOCAL node of green triangle, i.e. node against splitted edge
		// gpl2 is the second (counterclockwise) and gpl3 is the third
		// edit neighbor of triangle we are about to add
		gtn = _triangles[t].neighbors(gpl2);
		if (gtn != -1) // well, if it exists
			_triangles[gtn].neighbors(_triangles[gtn].neighbor2node(t)) = _triangles.size();
		_triangles.push_back(Triangle(
			_triangles[t].nodes(gpl1), splittingNode, _triangles[t].nodes(gpl3), // nodes counterclockwise
			first, gtn, t // neighbors
			));
		// ok, lets edit t
		_triangles[t]
			.nodes(_triangles[t].nodes(gpl1), _triangles[t].nodes(gpl2), splittingNode)
			.neighbors(second, _triangles.size() - 1, _triangles[t].neighbors(gpl3));
	};

	for (iter1 = indicies.begin(); iter1 != indicies.end(); ++iter1) {
		i = *iter1;
		// we will need ith nodes and neighbors later
		old = _triangles[i];
		p1 = old.nodes(0); p2 = old.nodes(1); p3 = old.nodes(2);
		t1 = old.neighbors(0); t2 = old.neighbors(1); t3 = old.neighbors(2);
		// add new points
		_nodes.push_back(_nodes[p1].midPoint(_nodes[p2]));
		_nodes.push_back(_nodes[p2].midPoint(_nodes[p3]));
		_nodes.push_back(_nodes[p3].midPoint(_nodes[p1]));
		// rp means red point, i.e. point added after red refinement
		rp1 = (rp2 = (rp3 = _nodes.size() - 1) - 1) - 1;
		// and we have also 3 red triangles yet to be added
		rt3 = (rt2 = (rt1 = _triangles.size()) + 1) + 1;
		// our ith triangle splits into 4 new ones
		// central one will take place of the old one
		_triangles[i]
			.nodes(rp1, rp2, rp3)
			// and its neighbors are known (3 other triangles) and will be added soon
			.neighbors(rt3, rt1, rt2);
		// now we have to add 3 other triangles
		// you have to draw them not to get confused w/ numeration
		// or just TRUST ME I AM A DOCTOR
		_triangles.push_back(Triangle(p1, rp1, rp3, i, -1, -1));
		_triangles.push_back(Triangle(rp1, p2, rp2, -1, i, -1));
		_triangles.push_back(Triangle(rp3, rp2, p3, -1, -1, i));
		// yet there 6 neighbors to be found
		// if we were dealing w/ bndry, we are gold (look at -1’s at previous 4 strings of code)
		// otherwise…
		if (t1 != -1) {
			iter2 = find(iter1, indicies.end(), t1);
			if (iter2 == indicies.end()) { // if we were not going to red-refine t1…
				// …then we will refine it green!
				// but first let’s deal w/ neighbors
				_triangles[rt2].neighbors(0) = _triangles.size(); // this one will be added soon
				_triangles[rt3].neighbors(0) = t1; // and t1 will be updated to fit in
				green(t1, rt2, rp2, rt3);
			}
			else {} // TODO
		}
		if (t2 != -1) {
			iter2 = find(iter1, indicies.end(), t2);
			if (iter2 == indicies.end()) {
				_triangles[rt3].neighbors(1) = _triangles.size();
				_triangles[rt1].neighbors(1) = t2;
				green(t2, rt3, rp3, rt1);
			}
			else {}
		}
		if (t3 != -1) {
			iter2 = find(iter1, indicies.end(), t3);
			if (iter2 == indicies.end()) {
				_triangles[rt1].neighbors(2) = _triangles.size();
				_triangles[rt2].neighbors(2) = t3;
				green(t3, rt1, rp1, rt2);
			}
			else {}
		}
	}
	return *this;
}