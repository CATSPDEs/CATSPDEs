#include "SymmetricContainer.hpp"
#include "Triangulation.hpp"
#include "array.hpp"
#include "SymmetricSkylineMatrix.hpp"
#include <set>

// our model problem:
//
// �nabla . (a nabla u) + cu = f,             if (x, y) in Omega,
// �n . (a nabla u) = kappa (u � g_D) � g_N,  if (x, y) in bndry of Omega,
//
// where a > 0, c >= c_0 > 0, kappa > 0, f, g_D, and g_N are given R � R �> R functions
//
// we convert our problem into a discrete one:
//
// (massMatrix + stiffnessMatrix + robinMatrix) . xi = loadVector + robinVector 
// (or A.xi = b for short),
//
// where matrix / element entries are given:
//
// massMatrix(i, j)      = dintt_Omega { c hatFunction_i hatFunction_j },
// stiffnessMatrix(i, j) = dintt_Omega { a nabla hatFunction_i . nabla hatFunction_j },
// robinMatrix(i, j)     = dintt_{bndry of Omega} { kappa hatFunction_i hatFunction_j },
// loadVector(i)         = dintt_Omega { f hatFunction_i },
// robinVector(i)        = dintt_{bndry of Omega} { (kappa g_D + g_N) hatFunction_i },
//
// where hatFunction_i denotes linear basis function taking unity on ith node and zero elsewhere
// so we have to solve n � n linear system, n := numb of nodes of the mesh Omega

inline double a(Node const& p) { // a > 0
	return 1.;
}

inline double c(Node const& p) { // c > c_0 >= 0
	return 1.;
}

inline double f(Node const& p) {
	return 1.;
}

inline double g_D(Node const& p) {
	return 0.;
}

inline double g_N(Node const& p) {
	return 0.;
}

inline double kappa(Node const& p) { // kappa > 0
	return 1.;
}

AdjacencyList generateAdjList(Triangulation const& Omega) {
	AdjacencyList adjList(Omega.numbOfNodes());
	for (size_t i = 0;i < Omega.numbOfTriangles();i++) {
		auto elementNodesIndicies = Omega.getNodesIndicies(i);
		sort(elementNodesIndicies.begin(), elementNodesIndicies.end());
		adjList[elementNodesIndicies[2]].insert(elementNodesIndicies[1]);
		adjList[elementNodesIndicies[2]].insert(elementNodesIndicies[0]);
		adjList[elementNodesIndicies[1]].insert(elementNodesIndicies[0]);
		//auto maxInd = *max_element(elementNodesIndicies.begin(), elementNodesIndicies.end());
		//for (auto ind : elementNodesIndicies)
		//	if (ind != maxInd)
		//		adjList[maxInd].insert(ind);
	}
	return adjList;
}

int main() {
	try {
		Triangulation Omega(Node(-1, -1), Node(1, 1), .99); // simple square mesh
		// data structures for final linear system A.xi = b
		auto adjList = generateAdjList(Omega);
		SymmetricSkylineMatrix A(adjList);
		vector<double> b(Omega.numbOfNodes(), 0), // load vector
		               xi(Omega.numbOfNodes()); // discrete solution			   
		SymmetricContainer<double> massMatrixLoc(3), // for hat functions on triangles 
		                           stiffnessMatrixLoc(3), // we have 3 � 3 element matricies
		                           robinMatrixLoc(2); // and 2 � 2 element matricies for Robin BCs (just like element matrix in 1D)
		array<double, 3> loadVectorLoc; // and their
		array<double, 2> robinVectorLoc; // friends, element vectors
		array<Node, 3> elementNodes; // nodes of current triangle
		Node leftNode, rightNode; // dummy nodes
		double elementArea, // area of ith triangle
		       edgeLength; // length of bndry edge of ith thiangle
		array<size_t, 3> l2g; // local to global mapping of nodes
		// dummy indicies
		size_t i;
		localIndex j, k;
		for (i = 0; i < Omega.numbOfTriangles(); ++i) {
			//
			// (1) quadratures over elements
			//
			// in order to assemble stiffness matrix and load vector,
			// it is convinient to iterate over mesh elements (i.e. triangles)
			elementNodes = Omega.getNodes(i); // get nodes of ith triangle
			elementArea = Omega.area(i); // compute area of ith triangle
			l2g = Omega.l2g(i); // local to global mapping of nodes of ith element
			// (1.1) compute local mass matrix
			massMatrixLoc(0, 0) = elementArea * (6 * c(elementNodes[0]) + 2 * c(elementNodes[1]) + 2 * c(elementNodes[2])) / 60.;
			massMatrixLoc(0, 1) = elementArea * (2 * c(elementNodes[0]) + 2 * c(elementNodes[1]) +     c(elementNodes[2])) / 60.;
			massMatrixLoc(0, 2) = elementArea * (3 * c(elementNodes[0]) +     c(elementNodes[1]) + 2 * c(elementNodes[2])) / 60.;
			massMatrixLoc(1, 1) = elementArea * (2 * c(elementNodes[0]) + 6 * c(elementNodes[1]) + 2 * c(elementNodes[2])) / 60.;
			massMatrixLoc(1, 2) = elementArea * (    c(elementNodes[0]) + 2 * c(elementNodes[1]) + 2 * c(elementNodes[2])) / 60.;
			massMatrixLoc(2, 2) = elementArea * (2 * c(elementNodes[0]) + 2 * c(elementNodes[1]) + 6 * c(elementNodes[2])) / 60.;
			// (1.2) compute local stiffness matrix
			stiffnessMatrixLoc(0, 0) = (elementNodes[1].x() - elementNodes[2].x()) * (elementNodes[1].x() - elementNodes[2].x()) + 
		                               (elementNodes[1].y() - elementNodes[2].y()) * (elementNodes[1].y() - elementNodes[2].y());
			stiffnessMatrixLoc(0, 1) = (elementNodes[0].x() - elementNodes[2].x()) * (elementNodes[2].x() - elementNodes[1].x()) +
			                           (elementNodes[0].y() - elementNodes[2].y()) * (elementNodes[2].y() - elementNodes[1].y());
			stiffnessMatrixLoc(0, 2) = (elementNodes[0].x() - elementNodes[1].x()) * (elementNodes[1].x() - elementNodes[2].x()) +
			                           (elementNodes[0].y() - elementNodes[1].y()) * (elementNodes[1].y() - elementNodes[2].y());
			stiffnessMatrixLoc(1, 1) = (elementNodes[0].x() - elementNodes[2].x()) * (elementNodes[0].x() - elementNodes[2].x()) +
			                           (elementNodes[0].y() - elementNodes[2].y()) * (elementNodes[0].y() - elementNodes[2].y());
			stiffnessMatrixLoc(1, 2) = (elementNodes[1].x() - elementNodes[0].x()) * (elementNodes[0].x() - elementNodes[2].x()) +
			                           (elementNodes[1].y() - elementNodes[0].y()) * (elementNodes[0].y() - elementNodes[2].y());
			stiffnessMatrixLoc(2, 2) = (elementNodes[0].x() - elementNodes[1].x()) * (elementNodes[0].x() - elementNodes[1].x()) +
			                           (elementNodes[0].y() - elementNodes[1].y()) * (elementNodes[0].y() - elementNodes[1].y());
			for (j = 0; j < 3; ++j)
				for (k = j; k < 3; ++k)
					stiffnessMatrixLoc(j, k) *= (a(elementNodes[0]) + a(elementNodes[1]) + a(elementNodes[2])) / elementArea / 12.;
			// (1.3) compute local load vector
			(loadVectorLoc = { 
				2 * f(elementNodes[0]) +     f(elementNodes[1]) + 2 * f(elementNodes[2]),
				    f(elementNodes[0]) + 2 * f(elementNodes[1]) +     f(elementNodes[2]),
				    f(elementNodes[0]) +     f(elementNodes[1]) + 2 * f(elementNodes[2])
			}) *= elementArea / 12.;
			// (1.4) assemble contributions
			for (j = 0; j < 3; ++j) {
				for (k = 0; k < 3; ++k)
					A(l2g[j], l2g[k]) += massMatrixLoc(j, k) + stiffnessMatrixLoc(j, k);
				b[l2g[j]] += loadVectorLoc[j];
			}
			//
			// (2) quadratures over edges
			//
			// iterate over list of local indicies of boundary nodes
			for (localIndex edgeIndex : Omega.getBoundaryIndicies(i)) {
				// if edgeIndex = 2, then the edge against second node of ith triangle
				// is part of the boundary
				// so we need to assemble BCs here
				edgeLength = Omega.length(i, edgeIndex);
				leftNode = elementNodes[nextIndex(edgeIndex)]; // these are nodes that 
				rightNode = elementNodes[nextIndex(nextIndex(edgeIndex))]; // define the edge
				// (2.1) compute local Robin matrix
				robinMatrixLoc(0, 0) = edgeLength * (3 * kappa(leftNode) +     kappa(rightNode)) / 12.;
				robinMatrixLoc(0, 1) = edgeLength * (    kappa(leftNode) +     kappa(rightNode)) / 12.;
				robinMatrixLoc(1, 1) = edgeLength * (    kappa(leftNode) + 3 * kappa(rightNode)) / 12.;
				// (2.2) compute local Robin vector
				(robinVectorLoc = {
					4 * g_N(leftNode) + 2 * g_N(rightNode) + 
					g_D(leftNode) * (3 * kappa(leftNode) + kappa(rightNode)) +
					g_D(rightNode) * (kappa(leftNode) + kappa(rightNode)),
					2 * g_N(leftNode) + 4 * g_N(rightNode) +
					g_D(leftNode) * (kappa(leftNode) + kappa(rightNode)) +
					g_D(rightNode) * (kappa(leftNode) + 3 * kappa(rightNode)),
				}) *= edgeLength / 12.;
				// (2.3) assemble contributions
				for (j = 0; j < 2; ++j) {
					for (k = 0; k < 2; ++k)
						A(l2g[j], l2g[k]) += robinMatrixLoc(j, k);
					b[l2g[j]] += robinVectorLoc[j];
				}
			}
		}
		xi = b / A; // compute our desrete solution
		A.save();
		cout << b << "\n\n";
		cout << xi << endl;
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}