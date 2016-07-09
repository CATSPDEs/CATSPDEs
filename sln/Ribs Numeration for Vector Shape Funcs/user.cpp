#include <fstream>
#include <string>
#include "Triangulation.hpp"

int main() {
	ifstream iNodes    ("Mathematica/Generate Mesh/nodes.dat"),
		     iTriangles("Mathematica/Generate Mesh/triangles.dat"),
		     iNeighbors("Mathematica/Generate Mesh/neighbors.dat");
	ofstream oNodes    ("Mathematica/Draw Mesh/nodes.dat"),
	         oTriangles("Mathematica/Draw Mesh/triangles.dat"),
	         oRibs     ("Mathematica/Draw Mesh/ribs.dat");
	try {
		Triangulation Omega(iNodes, iTriangles, iNeighbors);
		RibsNumeration ribs = Omega.computeRibsNumeration();
		Omega.save(oNodes, oTriangles);
		oRibs << ribs;
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}