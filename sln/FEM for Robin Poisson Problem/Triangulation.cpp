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