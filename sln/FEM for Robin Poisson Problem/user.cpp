#include <fstream>
#include "Triangulation.hpp"

int main() {
	ofstream redIndicies("Mathematica/tRed.dat");
	Triangulation K(Node(-2.5, -1), Node(2.5, 1), .2);
	cout << K.area(0) << ' ' << K.length(0, 0) << '\n'; 
	K.save(ofstream("Mathematica/nOld.dat"), ofstream("Mathematica/tOld.dat"));
	Indicies L;
	for (size_t i = 0; i < K.numbOfTriangles(); ++i)
		if (rand() % 2 && rand() % 2) {
			L.push_back(i); // ~25% of triangles will be red-refined
			redIndicies << i << '\n';
		}
	K.refine(L);
	K.save(ofstream("Mathematica/nNew.dat"), ofstream("Mathematica/tNew.dat"));
	// area of the very 1st triangle and length of its 1st side
	// (i.e. side against its 1st node)
	// it happens to be hypotenuse (and legs measure is equal to unity)
	return 0;
}