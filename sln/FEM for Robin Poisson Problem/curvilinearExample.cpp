#include <fstream>
#include "Triangulation.hpp"

inline Node parabolaCurve(double theta) {
	double s = 1 - 2 * theta;
	return Node(s, 1 - s * s);
}

double const PI = 3.141592653589793;

inline Node semiCircleCurve(double theta) { // unit speed parametrization
	double s = sin(PI / 2 * theta);
	s *= s;
	return Node(s - 1, -sqrt(.25 - (s - .5) * (s - .5)));
}

int main() {
	vector<Node> nodes = { Node(-1, 0), Node(0, 0), Node(0, 1), Node(1, 0), Node(-.5, -.5) };
	vector<Curve> curves = { parabolaCurve, semiCircleCurve };
	vector<CurvilinearEdge> edges = { 
		CurvilinearEdge(.5, 1., 0), 
		CurvilinearEdge(0., .5, 1),
		CurvilinearEdge(.5, 1., 1),
		CurvilinearEdge(0., .5, 0) 
	};
	vector<Triangle> triangles = {
		Triangle(0, 1, 2, 1, -2, 2),
		Triangle(1, 3, 2, -5, 0, -1),
		Triangle(0, 4, 1, -4, 0, -3)
	};
	Triangulation K(nodes, triangles, curves, edges);
	Indicies L = { 0, 1, 2 };
	K.refine(L);
	for (unsigned i = 0; i < 1; ++i) {
		while (L.back() < K.numbOfTriangles() - 1) L.push_back(L.back() + 1);
		K.refine(L);
	}
	K.save(ofstream("Mathematica/nCurve.dat"), ofstream("Mathematica/tCurve.dat"));
	return 0;
}