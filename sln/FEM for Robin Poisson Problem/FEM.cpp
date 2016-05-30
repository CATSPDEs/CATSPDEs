#include <fstream>
#include "SymmetricContainer.hpp" // for local matrices
#include "SymmetricCSlRMatrix.hpp" // for final linear system matrix
#include "krylov.hpp" // conjugate gradients
#include "array.hpp" // utility for array operations

// our model problem:
//
// –nabla . (a nabla u) + cu = f,             if (x, y) in Omega,
// –n . (a nabla u) = kappa (u – g_D) – g_N,  if (x, y) in bndry of Omega,
//
// where a > 0, c >= c_0 > 0, kappa > 0, f, g_D, and g_N are given R × R —> R functions
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
// so we have to solve n × n linear system, n := numb of nodes of the mesh Omega

// input R × R —> R functions (PDE and BCs) and the mesh (domain)
#include "PureDirichetProblem.hpp"

int main() {
	try {
		// data structures for final linear system A.xi = b:
		SymmetricCSlRMatrix A(Omega.generateAdjList()); // build final matrix portrait
		vector<double> b(Omega.numbOfNodes(), 0.), // load vector
		               xi(Omega.numbOfNodes(), 0.); // discrete solution	
		// data structures for assemby of A and b:
		SymmetricContainer<double> massMatrixLoc(3), // for hat functions on triangles 
		                           stiffnessMatrixLoc(3), // we have 3 × 3 element matricies
		                           robinMatrixLoc(2); // and 2 × 2 element matricies for Robin BCs (just like element matrix in 1D)
		array<double, 3> loadVectorLoc; // and their
		array<double, 2> robinVectorLoc; // friends, element vectors
		array<Node, 3> elementNodes, // nodes of current triangle
                       elementMiddleNodes; // and nodes on the middle of edges
		Node leftNode, rightNode; // dummy nodes
		double measure; // area of ith triangle / length of bndry edge of ith thiangle
		array<size_t, 3> l2g_elem; // local to global mapping of nodes on the element
		array<size_t, 2> l2g_edge; // and on the edge
		// dummy indicies
		size_t i;
		localIndex j, k, leftNodeIndex, rightNodeIndex;
		for (i = 0; i < Omega.numbOfTriangles(); ++i) {
			//
			// (1) quadratures over elements
			//
			// in order to assemble stiffness matrix and load vector,
			// it is convinient to iterate over mesh elements (i.e. triangles)
			elementNodes = Omega.getNodes(i); // get nodes of ith triangle
			for (j = 0; j < 3; ++j) // and middle nodes of its edges
				elementMiddleNodes[j] = elementNodes[k = nextIndex(j)].midPoint(elementNodes[nextIndex(k)]);
			measure = Omega.area(i); // compute area of ith triangle
			l2g_elem = Omega.l2g(i); // local to global mapping of nodes of ith element
			// (1.1) compute local mass matrix
			massMatrixLoc(0, 0) = measure * (6. * c(elementNodes[0]) + 2. * c(elementNodes[1]) + 2. * c(elementNodes[2])) / 60.;
			massMatrixLoc(0, 1) = measure * (2. * c(elementNodes[0]) + 2. * c(elementNodes[1]) +      c(elementNodes[2])) / 60.;
			massMatrixLoc(0, 2) = measure * (2. * c(elementNodes[0]) +      c(elementNodes[1]) + 2. * c(elementNodes[2])) / 60.;
			massMatrixLoc(1, 1) = measure * (2. * c(elementNodes[0]) + 6. * c(elementNodes[1]) + 2. * c(elementNodes[2])) / 60.;
			massMatrixLoc(1, 2) = measure * (     c(elementNodes[0]) + 2. * c(elementNodes[1]) + 2. * c(elementNodes[2])) / 60.;
			massMatrixLoc(2, 2) = measure * (2. * c(elementNodes[0]) + 2. * c(elementNodes[1]) + 6. * c(elementNodes[2])) / 60.;
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
					// quadratures calculated assuming a(x, y) lives in P_1(ith triangle), 
					// i.e. a(x, y) is linear combination of {x, y, 1}:
					// stiffnessMatrixLoc(j, k) *= (a(elementNodes[0]) + a(elementNodes[1]) + a(elementNodes[2])) / measure / 12.;
					// but we can do better:
					// quadratures calculated assuming a(x, y) lives in P_2(ith triangle), 
					// i.e. a(x, y) is linear combination of {x^2, y^2, xy, x, y, 1}:
					stiffnessMatrixLoc(j, k) *= (a(elementMiddleNodes[0]) + a(elementMiddleNodes[1]) + a(elementMiddleNodes[2])) / measure / 12.;
			// (1.3) compute local load vector
			(loadVectorLoc = { 
				2. * f(elementNodes[0]) +      f(elementNodes[1]) +      f(elementNodes[2]),
				     f(elementNodes[0]) + 2. * f(elementNodes[1]) +      f(elementNodes[2]),
				     f(elementNodes[0]) +      f(elementNodes[1]) + 2. * f(elementNodes[2])
			}) *= measure / 12.;
			// (1.4) assemble contributions
			for (j = 0; j < 3; ++j) {
				for (k = j; k < 3; ++k)
					A(l2g_elem[j], l2g_elem[k]) += massMatrixLoc(j, k) + stiffnessMatrixLoc(j, k);
				b[l2g_elem[j]] += loadVectorLoc[j];
			}
			//
			// (2) quadratures over edges
			//
			// iterate over list of local indicies of boundary nodes
			for (localIndex edgeIndex : Omega.getBoundaryIndicies(i)) {
				// if edgeIndex = 2, then the edge against second node of ith triangle
				// is part of the boundary
				// so we need to assemble BCs here
				leftNodeIndex = nextIndex(edgeIndex); // local indicies of nodes that
				rightNodeIndex = nextIndex(leftNodeIndex); // define the edge
				leftNode = elementNodes[leftNodeIndex]; // and the nodes 
				rightNode = elementNodes[rightNodeIndex]; // themselves
				l2g_edge[0] = l2g_elem[leftNodeIndex]; // local to global nodes
				l2g_edge[1] = l2g_elem[rightNodeIndex]; // numeration mapping 
				measure = Omega.length(i, edgeIndex);
				// (2.1) compute local Robin matrix
				robinMatrixLoc(0, 0) = measure * (3. * kappa(leftNode) +      kappa(rightNode)) / 12.;
				robinMatrixLoc(0, 1) = measure * (     kappa(leftNode) +      kappa(rightNode)) / 12.;
				robinMatrixLoc(1, 1) = measure * (     kappa(leftNode) + 3. * kappa(rightNode)) / 12.;
				// (2.2) compute local Robin vector
				(robinVectorLoc = {
					4. * g_N(leftNode) + 2. * g_N(rightNode) + 
					g_D(leftNode) * (3. * kappa(leftNode) + kappa(rightNode)) +
					g_D(rightNode) * (kappa(leftNode) + kappa(rightNode)),
					2. * g_N(leftNode) + 4. * g_N(rightNode) +
					g_D(leftNode) * (kappa(leftNode) + kappa(rightNode)) +
					g_D(rightNode) * (kappa(leftNode) + 3. * kappa(rightNode))
				}) *= measure / 12.;
				// (2.3) assemble contributions
				for (j = 0; j < 2; ++j) {
					for (k = j; k < 2; ++k)
						A(l2g_edge[j], l2g_edge[k]) += robinMatrixLoc(j, k);
					b[l2g_edge[j]] += robinVectorLoc[j];
				}
			}
		}
		// now we are ready to compute xi, A.xi = b
		xi = CG(A, b, xi, 10e-50);
		ofstream xiOutput("Mathematica/Model Problem Analysis/xi.dat");
		xiOutput << xi;
		Omega.save(ofstream("Mathematica/Model Problem Analysis/n.dat"), ofstream("Mathematica/Model Problem Analysis/t.dat"));
		
		vector<double> u(Omega.numbOfNodes());
		for (i = 0; i < u.size(); ++i)
			u[i] = 2. * Omega.getNode(i).x() + 3. * Omega.getNode(i).y() + 1.;
		ofstream bOutput("Mathematica/Model Problem Analysis/RHS.dat");
		bOutput.precision(15); // double precision
		bOutput << scientific << showpos;
		for (i = 0; i < u.size(); ++i)
			bOutput << b[i] << ' ' << (A * u)[i] << '\n';
	}
	catch (exception const & e) {
		cout << e.what() << endl;
	}
}