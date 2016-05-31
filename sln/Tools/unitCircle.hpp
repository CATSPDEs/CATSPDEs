#pragma once
#include "Triangulation.hpp"

// mesh for unit circle w/ center at origo

double const PI = 3.141592653589793;

inline Node circleCurve(double theta) { // unit speed parametrization
	return Node(cos(2 * PI * theta), sin(2 * PI * theta));
}

vector<Node> nodes = { Node(-1., 0.), Node(0., 0.), Node(0., 1.), Node(1., 0.), Node(0., -1.) };

vector<Curve> curves = { circleCurve };

vector<CurvilinearEdge> edges = {
	CurvilinearEdge(0., .25, 0),
	CurvilinearEdge(.25, .5, 0),
	CurvilinearEdge(.5, .75, 0),
	CurvilinearEdge(.75, 1., 0)
};

vector<Triangle> triangles = {
	Triangle(0, 1, 2, 3, -3, 1),
	Triangle(0, 4, 1, 2, 0, -4),
	Triangle(4, 3, 1, 3, 1, -5),
	Triangle(3, 2, 1, 0, 2, -2)
};

Triangulation Omega(nodes, triangles, curves, edges);
