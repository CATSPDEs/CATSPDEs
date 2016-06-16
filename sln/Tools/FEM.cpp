#include "FEM.hpp"
#include "SymmetricCSlRMatrix.hpp" // for final linear system matrix
#include "krylov.hpp" // conjugate gradients
#include "array.hpp" // utility for array operations

vector<double> FEM::computeDiscreteSolution(DiffusionReactionEqn const & PDE, 
                                            BoundaryConditions const & BCs, 
                                            Triangulation& Omega) {
	// data structures for final linear system A.xi = b:
	SymmetricCSlRMatrix A(Omega.generateAdjList()); // build final matrix portrait
	vector<double> b(Omega.numbOfNodes(), 0.), // load vector
	               xi(Omega.numbOfNodes(), 0.); // discrete solution	
	// data structures for assemby of A and b:
	SymmetricContainer<double> localMassMatrix(3), // for hat functions on triangles 
	                           localStiffnessMatrix(3), // we have 3 × 3 element matricies
	                           localRobinMatrix(2); // and 2 × 2 element matricies for Robin BCs (just like element matrix in 1D)
	array<double, 3> localLoadVector; // and their
	array<double, 2> localRobinVector; // friends, element vectors
	array<Node, 3> elementNodes, // nodes of the current triangle
				   elementMiddleNodes; // and nodes on the middle of edges
	array<Node, 2> edgeNodes; // nodes spanning an edge of the current triangle that is part of bndry
	double measure; // area of ith triangle / length of bndry edge of ith thiangle
	array<size_t, 3> l2g_elem; // local to global mapping of nodes on the element
	array<size_t, 2> l2g_edge; // and on the edge
	localIndex j, k, leftNodeIndex, rightNodeIndex; // dummy indicies
	for (size_t i = 0; i < Omega.numbOfTriangles(); ++i) {
		// (1) quadratures over elements
		// in order to assemble stiffness matrix and load vector,
		// it is convenient to iterate over mesh elements (i.e. triangles)
		elementNodes = Omega.getNodes(i); // get nodes of ith triangle
		for (j = 0; j < 3; ++j) // and middle nodes of its edges
			elementMiddleNodes[j] = elementNodes[k = nextIndex(j)].midPoint(elementNodes[nextIndex(k)]);
		measure = Omega.area(i); // compute area of ith triangle
		l2g_elem = Omega.l2g(i); // local to global mapping of nodes of ith element
		// compute
		// (1.1) local mass matrix,
		// (1.2) local stiffness matrix, and
		// (1.3) local load vector
		localStiffnessMatrix = computeLocalStiffnessMatrix(PDE.diffusionTerm(), elementNodes, elementMiddleNodes, measure);
		localMassMatrix      = computeLocalMassMatrix(PDE.reactionTerm(), elementNodes, measure);
		localLoadVector      = computeLocalLoadVector(PDE.forceTerm(), elementNodes, elementMiddleNodes, measure);
		// (1.4) assemble contributions
		for (j = 0; j < 3; ++j) {
			for (k = j; k < 3; ++k)
				A(l2g_elem[j], l2g_elem[k]) += localMassMatrix(j, k) + localStiffnessMatrix(j, k);
			b[l2g_elem[j]] += localLoadVector[j];
		}
		// (2) quadratures over edges
		// iterate over list of local indicies of boundary nodes
		for (localIndex edgeIndex : Omega.getBoundaryIndicies(i)) {
			// if edgeIndex = 2, then the edge against second node of ith triangle
			// is part of the boundary
			// so we need to assemble BCs here
			leftNodeIndex = nextIndex(edgeIndex); // local indicies of nodes that
			rightNodeIndex = nextIndex(leftNodeIndex); // define the edge
			edgeNodes = { elementNodes[leftNodeIndex], elementNodes[rightNodeIndex] }; // and the nodes themselves
			l2g_edge[0] = l2g_elem[leftNodeIndex]; // local to global nodes
			l2g_edge[1] = l2g_elem[rightNodeIndex]; // numeration mapping 
			measure = Omega.length(i, edgeIndex);
			// compute
			// (2.1) local Robin matrix
			// (2.2) local Robin vector
			localRobinMatrix = computeLocalRobinMatrix(BCs.RobinCoefficient(), edgeNodes, measure);
			localRobinVector = computeLocalRobinVector(BCs, edgeNodes, measure);
			// (2.3) assemble contributions
			for (j = 0; j < 2; ++j) {
				for (k = j; k < 2; ++k)
					A(l2g_edge[j], l2g_edge[k]) += localRobinMatrix(j, k);
				b[l2g_edge[j]] += localRobinVector[j];
			}
		}
	}
	// now we are ready to compute xi, A.xi = b
	xi = CG(A, b, xi, 10e-70);
	return xi;
}

SymmetricContainer<double> FEM::computeLocalMassMatrix(Function reactionTerm, 
	                                                   array<Node, 3>& nodes, 
	                                                   double area) {
	SymmetricContainer<double> m(3);
	m(0, 0) = area * (6. * reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
	m(0, 1) = area * (2. * reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + reactionTerm(nodes[2])) / 60.;
	m(0, 2) = area * (2. * reactionTerm(nodes[0]) + reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
	m(1, 1) = area * (2. * reactionTerm(nodes[0]) + 6. * reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
	m(1, 2) = area * (reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
	m(2, 2) = area * (2. * reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + 6. * reactionTerm(nodes[2])) / 60.;
	return m;
}

SymmetricContainer<double> FEM::computeLocalStiffnessMatrix(Function diffusionTerm,
                                                            array<Node, 3>& nodes,
													        array<Node, 3>& middleNodes,
                                                            double area) {
	SymmetricContainer<double> s(3);
	s(0, 0) = (nodes[1].x() - nodes[2].x()) * (nodes[1].x() - nodes[2].x()) +
	          (nodes[1].y() - nodes[2].y()) * (nodes[1].y() - nodes[2].y());
	s(0, 1) = (nodes[0].x() - nodes[2].x()) * (nodes[2].x() - nodes[1].x()) +
	          (nodes[0].y() - nodes[2].y()) * (nodes[2].y() - nodes[1].y());
	s(0, 2) = (nodes[0].x() - nodes[1].x()) * (nodes[1].x() - nodes[2].x()) +
	          (nodes[0].y() - nodes[1].y()) * (nodes[1].y() - nodes[2].y());
	s(1, 1) = (nodes[0].x() - nodes[2].x()) * (nodes[0].x() - nodes[2].x()) +
	          (nodes[0].y() - nodes[2].y()) * (nodes[0].y() - nodes[2].y());
	s(1, 2) = (nodes[1].x() - nodes[0].x()) * (nodes[0].x() - nodes[2].x()) +
	          (nodes[1].y() - nodes[0].y()) * (nodes[0].y() - nodes[2].y());
	s(2, 2) = (nodes[0].x() - nodes[1].x()) * (nodes[0].x() - nodes[1].x()) +
	          (nodes[0].y() - nodes[1].y()) * (nodes[0].y() - nodes[1].y());
	for (localIndex i = 0; i < 3; ++i)
		for (localIndex j = i; j < 3; ++j)
			// quadratures calculated assuming diffusionTerm(x, y) 
			// ( where (x, y) is a point from triangle spanned by nodes array )  
			// lives in P_1, 
			// i.e. diffusionTerm(x, y) is linear combination of {x, y, 1}:
			/*
			s(i, j) *= (diffusionTerm(nodes[0]) + diffusionTerm(nodes[1]) + diffusionTerm(nodes[2])) / area / 12.;
			*/
			// but we can do better:
			// quadratures calculated assuming diffusionTerm(x, y) lives in P_2(ith triangle), 
			// i.e. diffusionTerm(x, y) is linear combination of {x^2, y^2, xy, x, y, 1}:
			s(i, j) *= (diffusionTerm(middleNodes[0]) + diffusionTerm(middleNodes[1]) + diffusionTerm(middleNodes[2])) / area / 12.;
	return s;
}

array<double, 3> FEM::computeLocalLoadVector(Function forceTerm,
                                             array<Node, 3>& nodes,
                                             array<Node, 3>& middleNodes,
	                                         double area) {
	array<double, 3> f;
	// …assuming forceTerm(x, y) lives in P_1:
	/*
	return (f = {
	2. * forceTerm(nodes[0]) + forceTerm(nodes[1]) + forceTerm(nodes[2]),
	     forceTerm(nodes[0]) + 2. * forceTerm(nodes[1]) + forceTerm(nodes[2]),
	     forceTerm(nodes[0]) + forceTerm(nodes[1]) + 2. * forceTerm(nodes[2])
	}) *= area / 12.; 
	*/
	// …assuming forceTerm(x, y) lives in P_2(ith triangle):
	return (f = {
		2. * forceTerm(nodes[0]) - forceTerm(nodes[1]) - forceTerm(nodes[2]) +
		4. * forceTerm(middleNodes[0]) + 8. * (forceTerm(middleNodes[1]) + forceTerm(middleNodes[2])),
		2. * forceTerm(nodes[1]) - forceTerm(nodes[0]) - forceTerm(nodes[2]) +
		4. * (2. * (forceTerm(middleNodes[0]) + forceTerm(middleNodes[2])) + forceTerm(middleNodes[1])),
		2. * (forceTerm(nodes[2]) + 4. * (forceTerm(middleNodes[0]) + forceTerm(middleNodes[1])) + 2. * forceTerm(middleNodes[2])) -
		forceTerm(nodes[0]) - forceTerm(nodes[1])
	}) *= area / 60.;
}


SymmetricContainer<double> FEM::computeLocalRobinMatrix(Function RobinCoefficient,
                                                        array<Node, 2>& nodes,
	                                                    double length) {
	SymmetricContainer<double> r(2);
	r(0, 0) = length * ( 3. * RobinCoefficient(nodes[0]) +      RobinCoefficient(nodes[1]) ) / 12.;
	r(0, 1) = length * (      RobinCoefficient(nodes[0]) +      RobinCoefficient(nodes[1]) ) / 12.;
	r(1, 1) = length * (      RobinCoefficient(nodes[0]) + 3. * RobinCoefficient(nodes[1]) ) / 12.;
	return r;
}

array<double, 2> FEM::computeLocalRobinVector(BoundaryConditions const & BCs,
											  array<Node, 2>& nodes,
	                                          double length) {
	array<double, 2> r;
	return (r = {
		4. * BCs.NeumannValue(nodes[0]) + 2. * BCs.NeumannValue(nodes[1]) +
		BCs.DirichletCondition(nodes[0]) * ( 3. * BCs.RobinCoefficient(nodes[0]) + BCs.RobinCoefficient(nodes[1]) ) +
		BCs.DirichletCondition(nodes[1]) * ( BCs.RobinCoefficient(nodes[0]) + BCs.RobinCoefficient(nodes[1]) ),
		2. * BCs.NeumannValue(nodes[0]) + 4. * BCs.NeumannValue(nodes[1]) +
		BCs.DirichletCondition(nodes[0]) * ( BCs.RobinCoefficient(nodes[0]) + BCs.RobinCoefficient(nodes[1]) ) +
		BCs.DirichletCondition(nodes[1]) * ( BCs.RobinCoefficient(nodes[0]) + 3. * BCs.RobinCoefficient(nodes[1]) )
	}) *= length / 12.;
}