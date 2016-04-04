#include <fstream>
#include <ctime>
#include "Triangulation.hpp"

double diff(clock_t t1, clock_t t2) {
	return double(t2 - t1) / CLOCKS_PER_SEC;
}

int main() {
	clock_t begTime, endTime;
	ofstream redIndicies("Mathematica/tRed.dat");
	// (1) generating mesh 
	begTime = clock();
	Triangulation K(Node(-1, -1), Node(1, 1), .00448);
	endTime = clock();
	size_t e = K.numbOfTriangles(),
		   n = K.numbOfNodes();
	cout << "(1) Mesh generated . . . " << diff(begTime, endTime) << "s\n"
		 << "    Numb of elements before refinement: " << e << '\n'
		 << "    Numb of nodes / SLAE size before refinement: " << n << '\n';
	// (2) saving mesh
	begTime = clock();
	K.save(ofstream("Mathematica/nOld.dat"), ofstream("Mathematica/tOld.dat"));
	endTime = clock();
	cout << "(2) Mesh saved . . . " << diff(begTime, endTime) << "s\n";
	// (3) generating / saving red list
	Indicies L;
	begTime = clock();
	for (size_t i = 0; i < K.numbOfTriangles(); ++i) {
		if (rand() % 2) { // ~50% of triangles will be red-refined
			L.push_back(i);
			redIndicies << i << '\n';
		}
	}
	endTime = clock();
	size_t r = L.size();
	cout << "(3) Red list generated and saved . . . " << diff(begTime, endTime) << "s\n"
		 << "    Size of red list before refinement: " << r << " (" << 100 * r / double(e) << "% of " << e << ")\n";
	// (4) refining mesh
	begTime = clock();
	K.refine(L);
	endTime = clock();
	cout << "(4) Mesh refined . . . " << diff(begTime, endTime) << "s\n"
		 << "    Numb of elements after refinement: " << K.numbOfTriangles() << '\n'
		 << "    Numb of nodes / SLAE size after refinement: " << K.numbOfNodes() << '\n'
		 << "    Size of red list after refinement: " << L.size() << " (" << 100 * L.size() / double(e) << "% of " << e << ")\n";
	// (5) saving refined mesh
	begTime = clock();
	K.save(ofstream("Mathematica/nNew.dat"), ofstream("Mathematica/tNew.dat"));
	endTime = clock();
	cout << "(5) Refined mesh saved . . . " << diff(begTime, endTime) << "s\n\n";
	// to sum up
	cout << 100 * r / double(e) << "% of elements were to be refined, SLAE size changed from " << n << " to " << K.numbOfNodes() << "\n\n";
	return 0;
}