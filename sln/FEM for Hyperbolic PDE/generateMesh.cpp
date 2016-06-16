#pragma once
#include "Triangulation.hpp"

double const PI = 3.141592653589793;

inline Node circleCurve(double theta) { // unit speed parametrization
	return Node(cos(2 * PI * theta), -sin(2 * PI * theta));
}

Triangulation generateMesh() {
	Triangulation Omega(Node(-1.5, -1.5), Node(0., 1.5));
	vector<Node> nodes = { 
		Node(.5, -1.5), Node(1.5, -1.5),
		circleCurve(.125),
		Node(1.5, -.5),
		Node(1., 0.), Node(1.5, 0.),
		Node(1.5, .5),
		circleCurve(.875),
		Node(.5, 1.5), Node(1.5, 1.5) 
	};
	vector<Curve> curves = { circleCurve };
	vector<CurvilinearEdge> edges = {
		CurvilinearEdge(.0, .125, 0),
		CurvilinearEdge(.125, .25, 0),
		CurvilinearEdge(.75, .875, 0),
		CurvilinearEdge(.875, 1., 0)
	};
	vector<Triangle> triangles = {
		Triangle(3, 28, 7, 37, 5, -1),
		Triangle(28, 30, 7, -3, 36, 38),
		Triangle(28, 29, 30, 39, 37, -1),
		Triangle(29, 31, 30, 40, 38, -1),
		Triangle(30, 31, 32, 41, -2, 39),
		Triangle(31, 33, 32, 42, 40, -1),
		Triangle(32, 33, 34, -1, 43, 41),
		Triangle(32, 34, 35, 44, -5, 42),
		Triangle(34, 37, 35, 45, 43, -1),
		Triangle(35, 37, 36, -1, 46, 44),
		Triangle(35, 36, 23, 47, -4, 45),
		Triangle(23, 36, 27, -1, 35, 46)
	};
	Omega._nodes.insert(Omega._nodes.end(), nodes.begin(), nodes.end());
	Omega._triangles.insert(Omega._triangles.end(), triangles.begin(), triangles.end());
	Omega._curves.insert(Omega._curves.end(), curves.begin(), curves.end());
	Omega._edges.insert(Omega._edges.end(), edges.begin(), edges.end());
	Omega._makeNeighbors(5, 36);
	Omega._makeNeighbors(35, 47);
	Indicies L = { 37, 38, 39, 40, 43, 44, 45, 46 };
	Omega.refine(L);
	return Omega;
}