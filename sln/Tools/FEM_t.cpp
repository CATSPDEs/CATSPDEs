#include "FEM_t.hpp"
#include "SymmetricCSlRMatrix.hpp" // for final linear system matrix
#include "krylov.hpp" // conjugate gradients
#include "array.hpp" // utility for array operations

vector<vector<double>> FEM_t::CN3(HyperbolicPDE const & PDE,
	                             TimeFrames const & t,
	                             Triangulation& Omega,
	                             InitialConditions const & ICs,
	                             BoundaryConditions_t& BCs,
	                             Function_t u) {
	if (t[0] < 0.) throw invalid_argument("t0 < 0 — time cannot take negative values!");
	if (t.size() < 2) throw invalid_argument("CN3() needs at least 2 time frames to compute soln");
	// data structures for final linear system A.xi[m] = b:
	SymmetricCSlRMatrix A(Omega.generateAdjList()); // build final matrix portrait
	vector<double> b(Omega.numbOfNodes(), 0.); // load vector
	vector<vector<double>> xi(t.size(), b); // discrete solution—xi[m][i] is our solution at time = t[m] and at node_i
	// data structures for assemby of A and b:
	SymmetricContainer<double> localMassMatrixChi(3), localMassMatrixSigma(3), // for hat functions on triangles 
	                           localStiffnessMatrix(3), // we have 3 × 3 element matricies
	                           localRobinMatrix_m(2), // and 2 × 2 element matricies for Robin BCs (just like element matrix in 1D)
	                           localRobinMatrix_m_2(2); // computed on t[m] and t[m - 2]—for CN3–scheme
	array<double, 3> localLoadVector_m, localLoadVector_m_2; // and their
	array<double, 2> localRobinVector_m, localRobinVector_m_2; // friends, element vectors
	array<Node, 3> elementNodes, // nodes of the current triangle
	               elementMiddleNodes; // and nodes on the middle of edges
	array<Node, 2> edgeNodes; // nodes spanning an edge of the current triangle that is part of bndry
	Node midPoint;
	double measure; // area of ith triangle / length of bndry edge of ith thiangle
	array<size_t, 3> l2g_elem; // local to global mapping of nodes on the element
	array<size_t, 2> l2g_edge; // and on the edge
	localIndex j, k, leftNodeIndex, rightNodeIndex; // dummy indicies
	array<array<double, 2>, 3> eta; // time–basis
	// stencil:
	// t_(-2)-----t_(-1)---t_0
	// eta_i = unity at t_(-i), zero elsewhere, eta_i = a t^2 + b t + c
	// eta[i][j] means (j+1)’s derivative of eta_i computed at t_0
	// we will need these derivatives to approximate differential operators
	// (I) apply initial conditions to get xi[0] and xi[1]
	for (size_t i = 0; i < Omega.numbOfNodes(); ++i) 
		xi[0][i] = ICs.initialPosition(Omega.getNode(i));
	if (u == emptyFunc_t) for (size_t i = 0; i < Omega.numbOfNodes(); ++i)
		// normally, we have to use initial velocity to approximate xi[1]…
		xi[1][i] = xi[0][i] + ICs.initialVelocity(Omega.getNode(i)) * (t[1] - t[0]);
	else for (size_t i = 0; i < Omega.numbOfNodes(); ++i)
		// …but for model problems we have exact soln, so
		xi[1][i] = u(Omega.getNode(i), t[1]);
	// (II) solve for xi[m], m = 2, 3, … w/ Crank–Nicolson scheme
	for (size_t m = 2; m < t.size(); ++m) {
		// compute derivatives
		eta[2][0] = (t[m] - t[m - 1])                 / (t[m - 2] - t[m - 1]) / (t[m - 2] - t[m]);
		eta[2][1] = 2.                                / (t[m - 2] - t[m - 1]) / (t[m - 2] - t[m]);
		eta[1][0] = (t[m - 2] - t[m])                 / (t[m - 2] - t[m - 1]) / (t[m - 1] - t[m]);
		eta[1][1] = - 2.                              / (t[m - 2] - t[m - 1]) / (t[m - 1] - t[m]);
		eta[0][0] = (2. * t[m] - t[m - 2] - t[m - 1]) / (t[m - 2] - t[m]) / (t[m - 1] - t[m]);
		eta[0][1] = 2.                                / (t[m - 2] - t[m]) / (t[m - 1] - t[m]);
		// assemble SLDE
		for (size_t i = 0; i < Omega.numbOfTriangles(); ++i) {
			// (1) quadratures over elements
			// in order to assemble A and b,
			// it is convenient to iterate over mesh elements (i.e. triangles)
			elementNodes = Omega.getNodes(i); // get nodes of ith triangle
			for (j = 0; j < 3; ++j) // and middle nodes of its edges
				elementMiddleNodes[j] = elementNodes[k = nextIndex(j)].midPoint(elementNodes[nextIndex(k)]);
			measure = Omega.area(i); // compute area of ith triangle
			l2g_elem = Omega.l2g(i); // local to global mapping of nodes of ith element
			// local matrices and vectors
			localMassMatrixChi = computeLocalMassMatrix(PDE.chi(), elementNodes, measure);
			localMassMatrixSigma = computeLocalMassMatrix(PDE.sigma(), elementNodes, measure);
			localStiffnessMatrix = computeLocalStiffnessMatrix(PDE.diffusionTerm(), elementNodes, elementMiddleNodes, measure);
			localLoadVector_m = computeLocalLoadVector(PDE.forceTerm(), t[m], elementNodes, elementMiddleNodes, measure);
			localLoadVector_m_2 = computeLocalLoadVector(PDE.forceTerm(), t[m - 2], elementNodes, elementMiddleNodes, measure);
			// assemble contributions
			for (j = 0; j < 3; ++j) {
				for (k = 0; k < 3; ++k) {
					if (k >= j) // symmetric
						A(l2g_elem[j], l2g_elem[k]) +=
							eta[0][0] * localMassMatrixSigma(j, k) +
							eta[0][1] * localMassMatrixChi(j, k) +
							localStiffnessMatrix(j, k) / 2.;
					b[l2g_elem[j]] -= (
						localStiffnessMatrix(j, k) * xi[m - 2][l2g_elem[k]] / 2. +
						eta[1][0] * localMassMatrixSigma(j, k) * xi[m - 1][l2g_elem[k]] + // parabolic
						eta[2][0] * localMassMatrixSigma(j, k) * xi[m - 2][l2g_elem[k]] +
						eta[1][1] * localMassMatrixChi(j, k)   * xi[m - 1][l2g_elem[k]] + // hyperbolic
						eta[2][1] * localMassMatrixChi(j, k)   * xi[m - 2][l2g_elem[k]]
					);
				}
				b[l2g_elem[j]] += (localLoadVector_m[j] + localLoadVector_m_2[j]) / 2.;
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
				// define BCs to apply
				midPoint = edgeNodes[0].midPoint(edgeNodes[1]);
				BCs.defineBCsAt(midPoint);
				// compute
				// (2.1) local Robin matrix
				// (2.2) local Robin vector
				localRobinMatrix_m = FEM_t::computeLocalRobinMatrix(BCs, t[m], edgeNodes, measure);
				localRobinMatrix_m_2 = FEM_t::computeLocalRobinMatrix(BCs, t[m - 2], edgeNodes, measure);
				localRobinVector_m = FEM_t::computeLocalRobinVector(BCs, t[m], edgeNodes, measure);
				localRobinVector_m_2 = FEM_t::computeLocalRobinVector(BCs, t[m - 2], edgeNodes, measure);
				// (2.3) assemble contributions
				for (j = 0; j < 2; ++j) {
					for (k = 0; k < 2; ++k) {
						if (k >= j) // symmetric
							A(l2g_edge[j], l2g_edge[k]) += localRobinMatrix_m(j, k) / 2.;
						b[l2g_edge[j]] -= localRobinMatrix_m_2(j, k) * xi[m - 2][l2g_edge[k]] / 2.;
					}
					b[l2g_edge[j]] += (localRobinVector_m[j] + localRobinVector_m_2[j]) / 2.;
				}
			}
		}
		
		// now we are ready to compute xi[m], A.xi[m] = b
		//ofstream oA("Mathematica/Wave Interference Animation/m.dat");
		//A.save(oA);
		//return xi;

		xi[m] = CG(A, b, xi[m - 1], 10e-70);
		// clear A and b 
		A.setZero();
		fill(b.begin(), b.end(), 0.);
	}
	return xi;
}

vector<vector<double>> FEM_t::BDF3(HyperbolicPDE const & PDE,
	                             TimeFrames const & t,
	                             Triangulation& Omega,
	                             InitialConditions const & ICs,
	                             BoundaryConditions_t& BCs,
	                             Function_t u) {
	if (t.size() < 2) throw invalid_argument("BDF3() needs at least 2 time frames to compute soln");
	// data structures for final linear system A.xi[m] = b:
	SymmetricCSlRMatrix A(Omega.generateAdjList()); // build final matrix portrait
	vector<double> b(Omega.numbOfNodes(), 0.); // load vector
	vector<vector<double>> xi(t.size(), b); // discrete solution—xi[m][i] is our solution at time = t[m] and at node_i
	// data structures for assemby of A and b:
	SymmetricContainer<double> localMassMatrixChi(3), localMassMatrixSigma(3), // for hat functions on triangles 
		                       localStiffnessMatrix(3), // we have 3 × 3 element matricies
		                       localRobinMatrix(2); // and 2 × 2 element matricies for Robin BCs (just like element matrix in 1D)
	array<double, 3> localLoadVector; // and their
	array<double, 2> localRobinVector; // friends, element vectors
	array<Node, 3> elementNodes, // nodes of the current triangle
		           elementMiddleNodes; // and nodes on the middle of edges
	array<Node, 2> edgeNodes; // nodes spanning an edge of the current triangle that is part of bndry
	Node midPoint;
	double measure; // area of ith triangle / length of bndry edge of ith thiangle
	array<size_t, 3> l2g_elem; // local to global mapping of nodes on the element
	array<size_t, 2> l2g_edge; // and on the edge
	localIndex j, k, leftNodeIndex, rightNodeIndex; // dummy indicies
	array<array<double, 2>, 3> eta; // time–basis
	// stencil:
	// t_(-2)-----t_(-1)---t_0
	// eta_i = unity at t_(-i), zero elsewhere, eta_i = a t^2 + b t + c
	// eta[i][j] means (j+1)’s derivative of eta_i computed at t_0
	// we will need these derivatives to approximate differential operators
	// (I) apply initial conditions to get xi[0] and xi[1]
	for (size_t i = 0; i < Omega.numbOfNodes(); ++i)
		xi[0][i] = ICs.initialPosition(Omega.getNode(i));
	if (u == emptyFunc_t) for (size_t i = 0; i < Omega.numbOfNodes(); ++i)
		// normally, we have to use initial velocity to approximate xi[1]…
		xi[1][i] = xi[0][i] + ICs.initialVelocity(Omega.getNode(i)) * (t[1] - t[0]);
	else for (size_t i = 0; i < Omega.numbOfNodes(); ++i)
		// …but for model problems we have exact soln, so
		xi[1][i] = u(Omega.getNode(i), t[1]);
	// (II) solve for xi[m], m = 2, 3, … w/ Crank–Nicolson scheme
	for (size_t m = 2; m < t.size(); ++m) {
		// compute derivatives
		eta[2][0] = (t[m] - t[m - 1]) / (t[m - 2] - t[m - 1]) / (t[m - 2] - t[m]);
		eta[2][1] = 2. / (t[m - 2] - t[m - 1]) / (t[m - 2] - t[m]);
		eta[1][0] = (t[m - 2] - t[m]) / (t[m - 2] - t[m - 1]) / (t[m - 1] - t[m]);
		eta[1][1] = -2. / (t[m - 2] - t[m - 1]) / (t[m - 1] - t[m]);
		eta[0][0] = (2. * t[m] - t[m - 2] - t[m - 1]) / (t[m - 2] - t[m]) / (t[m - 1] - t[m]);
		eta[0][1] = 2. / (t[m - 2] - t[m]) / (t[m - 1] - t[m]);
		// assemble SLDE
		for (size_t i = 0; i < Omega.numbOfTriangles(); ++i) {
			// (1) quadratures over elements
			// in order to assemble A and b,
			// it is convenient to iterate over mesh elements (i.e. triangles)
			elementNodes = Omega.getNodes(i); // get nodes of ith triangle
			for (j = 0; j < 3; ++j) // and middle nodes of its edges
				elementMiddleNodes[j] = elementNodes[k = nextIndex(j)].midPoint(elementNodes[nextIndex(k)]);
			measure = Omega.area(i); // compute area of ith triangle
			l2g_elem = Omega.l2g(i); // local to global mapping of nodes of ith element
			// local matrices and vectors
			localMassMatrixChi = computeLocalMassMatrix(PDE.chi(), elementNodes, measure);
			localMassMatrixSigma = computeLocalMassMatrix(PDE.sigma(), elementNodes, measure);
			localStiffnessMatrix = computeLocalStiffnessMatrix(PDE.diffusionTerm(), elementNodes, elementMiddleNodes, measure);
			localLoadVector = computeLocalLoadVector(PDE.forceTerm(), t[m], elementNodes, elementMiddleNodes, measure);
			// assemble contributions
			for (j = 0; j < 3; ++j) {
				for (k = 0; k < 3; ++k) {
					if (k >= j) // symmetric
						A(l2g_elem[j], l2g_elem[k]) +=
							eta[0][0] * localMassMatrixSigma(j, k) +
							eta[0][1] * localMassMatrixChi(j, k) +
							localStiffnessMatrix(j, k);
					b[l2g_elem[j]] -= (
							eta[1][0] * localMassMatrixSigma(j, k) * xi[m - 1][l2g_elem[k]] + // parabolic
							eta[2][0] * localMassMatrixSigma(j, k) * xi[m - 2][l2g_elem[k]] +
							eta[1][1] * localMassMatrixChi(j, k)   * xi[m - 1][l2g_elem[k]] + // hyperbolic
							eta[2][1] * localMassMatrixChi(j, k)   * xi[m - 2][l2g_elem[k]]
						);
				}
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
				// define BCs to apply
				midPoint = edgeNodes[0].midPoint(edgeNodes[1]);
				BCs.defineBCsAt(midPoint);
				// compute
				// (2.1) local Robin matrix
				// (2.2) local Robin vector
				localRobinMatrix = FEM_t::computeLocalRobinMatrix(BCs, t[m], edgeNodes, measure);
				localRobinVector = FEM_t::computeLocalRobinVector(BCs, t[m], edgeNodes, measure);
				// (2.3) assemble contributions
				for (j = 0; j < 2; ++j) {
					for (k = j; k < 2; ++k) 
						A(l2g_edge[j], l2g_edge[k]) += localRobinMatrix(j, k);
					b[l2g_edge[j]] += localRobinVector[j];
				}
			}
		}
		// now we are ready to compute xi[m], A.xi[m] = b
		//ofstream oA("Mathematica/Model Problem Analysis/m.dat");
		//A.save(oA);
		xi[m] = CG(A, b, xi[m - 1], 10e-70);
		// clear A and b 
		A.setZero();
		fill(b.begin(), b.end(), 0.);
	}
	return xi;
}


SymmetricContainer<double> FEM_t::computeLocalMassMatrix(Function reactionTerm,
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

SymmetricContainer<double> FEM_t::computeLocalStiffnessMatrix(Function diffusionTerm,
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

array<double, 3> FEM_t::computeLocalLoadVector(Function_t forceTerm,
                                              double t,
                                              array<Node, 3>& nodes,
                                              array<Node, 3>& middleNodes,
                                              double area) {
	array<double, 3> f;
	return (f = {
		2. * forceTerm(nodes[0], t) - forceTerm(nodes[1], t) - forceTerm(nodes[2], t) +
		4. * forceTerm(middleNodes[0], t) + 8. * (forceTerm(middleNodes[1], t) + forceTerm(middleNodes[2], t)),
		2. * forceTerm(nodes[1], t) - forceTerm(nodes[0], t) - forceTerm(nodes[2], t) +
		4. * (2. * (forceTerm(middleNodes[0], t) + forceTerm(middleNodes[2], t)) + forceTerm(middleNodes[1], t)),
		2. * (forceTerm(nodes[2], t) + 4. * (forceTerm(middleNodes[0], t) + forceTerm(middleNodes[1], t)) + 2. * forceTerm(middleNodes[2], t)) -
		forceTerm(nodes[0], t) - forceTerm(nodes[1], t)
	}) *= area / 60.;
}


SymmetricContainer<double> FEM_t::computeLocalRobinMatrix(BoundaryConditions_t const & BCs,
                                                         double t,
                                                         array<Node, 2>& nodes,
	                                                     double length) {
	SymmetricContainer<double> r(2);
	r(0, 0) = length * (3. * BCs.RobinCoefficient(nodes[0], t) + BCs.RobinCoefficient(nodes[1], t)) / 12.;
	r(0, 1) = length * (BCs.RobinCoefficient(nodes[0], t) + BCs.RobinCoefficient(nodes[1], t)) / 12.;
	r(1, 1) = length * (BCs.RobinCoefficient(nodes[0], t) + 3. * BCs.RobinCoefficient(nodes[1], t)) / 12.;
	return r;
}

array<double, 2> FEM_t::computeLocalRobinVector(BoundaryConditions_t const & BCs,
                                               double t,
                                               array<Node, 2>& nodes,
                                               double length) {
	array<double, 2> r;
	return (r = {
		4. * BCs.NeumannValue(nodes[0], t) + 2. * BCs.NeumannValue(nodes[1], t) +
		BCs.DirichletCondition(nodes[0], t) * (3. * BCs.RobinCoefficient(nodes[0], t) + BCs.RobinCoefficient(nodes[1], t)) +
		BCs.DirichletCondition(nodes[1], t) * (BCs.RobinCoefficient(nodes[0], t) + BCs.RobinCoefficient(nodes[1], t)),
		2. * BCs.NeumannValue(nodes[0], t) + 4. * BCs.NeumannValue(nodes[1], t) +
		BCs.DirichletCondition(nodes[0], t) * (BCs.RobinCoefficient(nodes[0], t) + BCs.RobinCoefficient(nodes[1], t)) +
		BCs.DirichletCondition(nodes[1], t) * (BCs.RobinCoefficient(nodes[0], t) + 3. * BCs.RobinCoefficient(nodes[1], t))
	}) *= length / 12.;
}

vector<vector<double>> FEM_t::constructVector(Function_t u, TimeFrames const & time, Triangulation& Omega) {
	vector<vector<double>> uVec(time.size(), vector<double>(Omega.numbOfNodes()));
	for (size_t m = 0; m < time.size(); ++m)
		for (size_t i = 0; i < Omega.numbOfNodes(); ++i)
			uVec[m][i] = u(Omega.getNode(i), time[m]);
	return uVec;
}