#include <fstream>
#include <numeric> // accumulate()
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
	Triangulation M(K);
	vector<double> q;
	double mean, minElement;
	ofstream quality("Mathematica/quality.dat");
	// (1) “good” mesh example
	Indicies L = { 0, 1 }; // smart hack — we will not get “bad” triangles now
	K.refine(L);
	K.save(ofstream("Mathematica/nGoodIni.dat"), ofstream("Mathematica/tGoodIni.dat"));
	K.refine(3); // refine 3 times
	K.save(ofstream("Mathematica/nGood.dat"), ofstream("Mathematica/tGood.dat"));
	// quality measures
	q = K.qualityMeasure();
	minElement = *min_element(q.begin(), q.end());
	mean = accumulate(q.begin(), q.end(), 0.) / q.size();
	quality << minElement << '\n' << mean << '\n';
	// (2) “bad” example
	M.refine();
	M.save(ofstream("Mathematica/nBadIni.dat"), ofstream("Mathematica/tBadIni.dat"));
	M.refine(3); 
	M.save(ofstream("Mathematica/nBad.dat"), ofstream("Mathematica/tBad.dat"));
	q = M.qualityMeasure();
	minElement = *min_element(q.begin(), q.end());
	mean = accumulate(q.begin(), q.end(), 0.) / q.size();
	quality << minElement << '\n' << mean;
	return 0;
}