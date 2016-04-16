#include "DenseMatrix.hpp"
#include "SymmetricContainer.hpp"
#include "Triangulation.hpp"
#include "array.hpp"

inline double a(Node const& p) {
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

inline double kappa(Node const& p) {
	return 1.;
}

int main() {
	try {
		Triangulation Omega(Node(-1, -1), Node(1, 1), .99); // simple square mesh
		// stiffness matrix and load vector data structures
		DenseMatrix stiffnessMatrix(Omega.numbOfNodes()); // in dense format (for simplicity—will be changed soon)
		SymmetricContainer<double> stiffnessMatrixLoc(3); // for hat functions on triangles we have 3 × 3 element matricies
		SymmetricContainer<double> robinMatrixLoc(2); // and 2 × 2 element matricies for Robin BCs (just like element matrix in 1D)
		vector<double> loadVector(Omega.numbOfNodes(), 0); // load vector
		array<double, 3> loadVectorLoc; // and its
		array<double, 2> robinVectorLoc; // element friends
		array<double, 3> b, c; // let j be local index on ith element, j in {0, 1, 2},
		// then we have 3 non-zero hat functions defined on this triangle (jth function takes 1 on jth node):
		// hatFunction_j (x, y) = a_j + b_j x + c_j y,
		// so the gradient of jth function is
		// nabla . hatFunction_j = <b_j, c_j> —
		// so <b[0], c[0]> is gradient of 0th function (and so on)
		// we will need this gradients to assemble stiffness matrix BTW
		array<Node, 3> elementNodes; // nodes of current triangle
		Node leftNode, rightNode; // dummy nodes
		double elementArea; // area of ith triangle
		// dummy indicies
		size_t i;
		localIndex j, k;
		// and list of boundary indicies (empty for interior elementa)
		LocalIndicies boundaryLocalIndicies;
		for (i = 0; i < Omega.numbOfTriangles(); ++i) {
			// in order to assemble stiffness matrix and load vector,
			// it is convinient to iterate over mesh elements (i.e. triangles)
			// get nodes of ith triangle
			elementNodes = Omega.getNodes(i);
			// compute area of ith triangle
			elementArea = Omega.area(i);
			// compute gradient components of hat–functions on ith triangle
			(b = { (elementNodes[1] - elementNodes[2]).y(), (elementNodes[2] - elementNodes[0]).y(), (elementNodes[0] - elementNodes[1]).y() }) /= 2 * elementArea;
			(c = { (elementNodes[2] - elementNodes[1]).x(), (elementNodes[0] - elementNodes[2]).x(), (elementNodes[2] - elementNodes[1]).x() }) /= 2 * elementArea;
			// compute local load vector
			(loadVectorLoc = { f(elementNodes[0]), f(elementNodes[1]), f(elementNodes[2]) }) *= elementArea / 3;
			// compute local stiffness matrix
			for (j = 0; j < 3; ++j)
				for (k = j; k < 3; ++k)
					stiffnessMatrixLoc(j, k) = a((elementNodes[0] + elementNodes[1] + elementNodes[2]) / 3) * (b[j] * b[k] + c[j] * c[k]) * Omega.area(i);
			// iterate over list of local indicies of boundary nodes
			for (localIndex edgeIndex : Omega.getBoundaryIndicies(i)) {
				// if edgeIndex = 2, then the edge against second node of ith triangle
				// is part of the boundary
				// so we need to assemble BCs here
				leftNode = elementNodes[nextIndex(edgeIndex)];
				rightNode = elementNodes[nextIndex(nextIndex(edgeIndex))];
				// nodes that make the edge
				// assemble Robin local matrix
				robinMatrixLoc(0, 0) = Omega.length(i, edgeIndex) * (kappa(leftNode) + 3 * kappa(rightNode)) / 12;
				robinMatrixLoc(0, 1) = Omega.length(i, edgeIndex) * (kappa(leftNode) + kappa(rightNode)) / 12;
				robinMatrixLoc(1, 1) = Omega.length(i, edgeIndex) * (3 * kappa(leftNode) + kappa(rightNode)) / 12;
				// assemble Robin local vector
				robinVectorLoc = {
					Omega.length(i, edgeIndex) * (
						(kappa(leftNode) + kappa(rightNode)) * g_D(leftNode) +
						(kappa(leftNode) + 3 * kappa(rightNode)) * g_D(rightNode) +
						2 * (g_N(leftNode) + 2 * g_N(rightNode))
					),
					Omega.length(i, edgeIndex) * (
						(3 * kappa(leftNode) + kappa(rightNode)) * g_D(leftNode) +
						(kappa(leftNode) + kappa(rightNode)) * g_D(rightNode) +
						2 * (2 * g_N(leftNode) + g_N(rightNode))
					)
				};
			}
			// all the conributions are computed,
			// so we finally ready to assemble global stiffness matrix
			// and load vector
			for (j = 0; j < 3; ++j) {
				for (k = 0; k < 3; ++k)
					stiffnessMatrix(Omega.l2g(i)[j], Omega.l2g(i)[k]) += stiffnessMatrixLoc(j, k);
				loadVector[Omega.l2g(i)[j]] += loadVectorLoc[j];
			}
			for (j = 0; j < 2; ++j) {
				for (k = 0; k < 2; ++k) 
					stiffnessMatrix(Omega.l2g(i)[j], Omega.l2g(i)[k]) += robinMatrixLoc(j, k);
				loadVector[Omega.l2g(i)[j]] += robinVectorLoc[j];
			}
		}
		stiffnessMatrix.save();
		cout << loadVector << "\n\n";
		cout << loadVector / stiffnessMatrix << endl;
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}