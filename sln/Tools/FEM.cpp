#include "FEM.hpp"
#include "SymmetricMatrix.hpp" // for local matrices
#include "DenseSquareMatrix.hpp"

extern SingletonLogger& logger;

namespace FEM {

	namespace DivGrad {

		boost::tuple<
			SymmetricCSlCMatrix<double>, // system matrix
			std::vector<double> // rhs vector
		> FEM::DivGrad::linearLagrangeAssembler(
			DiffusionReactionEqn2D const & PDE, 
			Triangulation const & Omega,
			DirichletScalarCondition2D const & DirichletBC
		) {
			// data structures for final linear system A.xi = b:
			SymmetricCSlCMatrix<double> A; // system matrix
			logger.beg("generate matrix pattern");
				A.generatePatternFrom(Omega.generateAdjList()); // compute matrix pattern
				logger.log("numb of nnzeros = " + std::to_string(A.nnz()));
			logger.end();
			std::vector<double> b(Omega.numbOfNodes(), 0.), // load vector
								xi(Omega.numbOfNodes(), 0.); // discrete solution	
			Triangle elementNodes, elementMiddleNodes; // …of current element
			std::array<Node2D, 2> ribNodes; // nodes spanning a rib of the current triangle that is part of bndry
			double measure; // area of ith triangle / length of bndry edge of ith thiangle
			std::unordered_map<Index, double> ind2val; // for Dirichlet BCs

			auto computeLocalStiffnessMatrix = [&]() { // we have 3 × 3 element matricies for linear Lagrange shapes
				SymmetricMatrix<double> s(3);
				s(0, 0) = (elementNodes[1][0] - elementNodes[2][0]) * (elementNodes[1][0] - elementNodes[2][0]) +
						  (elementNodes[1][1] - elementNodes[2][1]) * (elementNodes[1][1] - elementNodes[2][1]);
				s(0, 1) = (elementNodes[0][0] - elementNodes[2][0]) * (elementNodes[2][0] - elementNodes[1][0]) +
						  (elementNodes[0][1] - elementNodes[2][1]) * (elementNodes[2][1] - elementNodes[1][1]);
				s(0, 2) = (elementNodes[0][0] - elementNodes[1][0]) * (elementNodes[1][0] - elementNodes[2][0]) +
						  (elementNodes[0][1] - elementNodes[1][1]) * (elementNodes[1][1] - elementNodes[2][1]);
				s(1, 1) = (elementNodes[0][0] - elementNodes[2][0]) * (elementNodes[0][0] - elementNodes[2][0]) +
						  (elementNodes[0][1] - elementNodes[2][1]) * (elementNodes[0][1] - elementNodes[2][1]);
				s(1, 2) = (elementNodes[1][0] - elementNodes[0][0]) * (elementNodes[0][0] - elementNodes[2][0]) +
						  (elementNodes[1][1] - elementNodes[0][1]) * (elementNodes[0][1] - elementNodes[2][1]);
				s(2, 2) = (elementNodes[0][0] - elementNodes[1][0]) * (elementNodes[0][0] - elementNodes[1][0]) +
						  (elementNodes[0][1] - elementNodes[1][1]) * (elementNodes[0][1] - elementNodes[1][1]);
				// quadratures calculated assuming diffusionTerm(x, y) lives in P_2(ith triangle), 
				// i.e. diffusionTerm(x, y) is linear combination of {x^2, y^2, xy, x, y, 1}:
				for (LocalIndex i = 0; i < 3; ++i)
					for (LocalIndex j = i; j < 3; ++j)
						s(i, j) *= (PDE.diffusionTerm(elementMiddleNodes[0]) + PDE.diffusionTerm(elementMiddleNodes[1]) + PDE.diffusionTerm(elementMiddleNodes[2])) / measure / 12.;
				return s;
			};
			auto computeLocalMassMatrix = [&]() {
				SymmetricMatrix<double> m(3);
				m(0, 0) = measure * (6. * PDE.reactionTerm(elementNodes[0]) + 2. * PDE.reactionTerm(elementNodes[1]) + 2. * PDE.reactionTerm(elementNodes[2])) / 60.;
				m(0, 1) = measure * (2. * PDE.reactionTerm(elementNodes[0]) + 2. * PDE.reactionTerm(elementNodes[1]) + PDE.reactionTerm(elementNodes[2])) / 60.;
				m(0, 2) = measure * (2. * PDE.reactionTerm(elementNodes[0]) + PDE.reactionTerm(elementNodes[1]) + 2. * PDE.reactionTerm(elementNodes[2])) / 60.;
				m(1, 1) = measure * (2. * PDE.reactionTerm(elementNodes[0]) + 6. * PDE.reactionTerm(elementNodes[1]) + 2. * PDE.reactionTerm(elementNodes[2])) / 60.;
				m(1, 2) = measure * (     PDE.reactionTerm(elementNodes[0]) + 2. * PDE.reactionTerm(elementNodes[1]) + 2. * PDE.reactionTerm(elementNodes[2])) / 60.;
				m(2, 2) = measure * (2. * PDE.reactionTerm(elementNodes[0]) + 2. * PDE.reactionTerm(elementNodes[1]) + 6. * PDE.reactionTerm(elementNodes[2])) / 60.;
				return m;
			};
			auto computeLocalLoadVector = [&]() {
				std::array<double, 3> f;
				// assuming forceTerm(x, y) lives in P_2(ith triangle):
				return (f = {
					2. * PDE.forceTerm(elementNodes[0]) - PDE.forceTerm(elementNodes[1]) - PDE.forceTerm(elementNodes[2]) +
					4. * PDE.forceTerm(elementMiddleNodes[0]) + 8. * (PDE.forceTerm(elementMiddleNodes[1]) + PDE.forceTerm(elementMiddleNodes[2])),
					2. * PDE.forceTerm(elementNodes[1]) - PDE.forceTerm(elementNodes[0]) - PDE.forceTerm(elementNodes[2]) +
					4. * (2. * (PDE.forceTerm(elementMiddleNodes[0]) + PDE.forceTerm(elementMiddleNodes[2])) + PDE.forceTerm(elementMiddleNodes[1])),
					2. * (PDE.forceTerm(elementNodes[2]) + 4. * (PDE.forceTerm(elementMiddleNodes[0]) + PDE.forceTerm(elementMiddleNodes[1])) + 2. * PDE.forceTerm(elementMiddleNodes[2])) -
					PDE.forceTerm(elementNodes[0]) - PDE.forceTerm(elementNodes[1])
				}) *= measure / 60.;
			};
			// we have 2 × 2 element matricies for Robin BCs (just like element matrix in 1D)
	
			logger.beg("assemble system matrix and rhs vector");
				for (size_t i = 0; i < Omega.numbOfElements(); ++i) {
					// (1) quadratures over elements
					// in order to assemble stiffness matrix and load vector,
					// it is convenient to iterate over mesh elements (i.e. triangles)
					elementNodes = Omega.getElement(i), // nodes of the current triangle
					elementMiddleNodes = midNodes(elementNodes); // and nodes on the middle of edges
					measure = area(elementNodes); // compute area of ith triangle
					// compute matrices and vectors
					auto localStiffnessMatrix = computeLocalStiffnessMatrix();
					auto localMassMatrix      = computeLocalMassMatrix();
					auto localLoadVector      = computeLocalLoadVector();
					// (1.4) assemble contributions
					auto l2g_elem = Omega.getNodesIndicies(i); // local to global mapping of nodes of ith element
					for (LocalIndex j : { 0, 1, 2 }) {
						for (LocalIndex k = j; k < 3; ++k)
							A(l2g_elem[j], l2g_elem[k]) += localMassMatrix(j, k) + localStiffnessMatrix(j, k);
						b[l2g_elem[j]] += localLoadVector[j];
					}
					// (2) quadratures over boundary edges
					for (LocalIndex j : { 0, 1, 2 }) {
						if (Omega.getNeighborsIndicies(i)[j] >= 0) continue;
						auto leftNodeLocalIndex = nextIndex(j), // local indicies of nodes that
							 rightNodeLocalIndex = nextIndex(leftNodeLocalIndex); // span the rib
						ribNodes = { elementNodes[leftNodeLocalIndex], elementNodes[rightNodeLocalIndex] }; // and the nodes themselves
						//measure = norm(edgeNodes[1] - edgeNodes[0]);
						//// define BCs to apply
						//midPoint = edgeNodes[0].midPoint(edgeNodes[1]);
						//BCs.defineBCsAt(midPoint);
						//// compute local matrix and vector
						//localRobinMatrix = computeLocalRobinMatrix(BCs, edgeNodes, measure);
						//localRobinVector = computeLocalRobinVector(BCs, edgeNodes, measure);
						//// (2.3) assemble contributions
						std::array<Index, 2> l2g_edge = { l2g_elem[leftNodeLocalIndex], l2g_elem[rightNodeLocalIndex] }; // local to global nodes numeration mapping 
						//for (j = 0; j < 2; ++j) {
						//	for (k = j; k < 2; ++k)
						//		A(l2g_edge[j], l2g_edge[k]) += localRobinMatrix(j, k);
						//	b[l2g_edge[j]] += localRobinVector[j];
						//}
						for (LocalIndex j : { 0, 1 })
							if (DirichletBC.shouldBeEnforcedAt(ribNodes[j])) 
								ind2val[l2g_edge[j]] = DirichletBC(ribNodes[j]);
					}
				}
			logger.end();
			logger.log("system size = " + std::to_string(A.getOrder()));
			logger.beg("enforce Dirichlet condition");
				A.enforceDirichletBCs(ind2val, b);
			logger.end();
			return boost::make_tuple(A, b);
		}

		namespace Multigrid {

			std::function<std::vector<double>(
				SymmetricCSlCMatrix<double>&, // system matrix
				std::vector<double> const &, // rhs
				std::vector<double> const & // initial guess
			)> smoother;
			Index gamma = 1; // numb of recursive coarse iterations (1 for V–cycle, 2 for W–cycle)
			std::vector<CSCMatrix<double>> prolongations; // prolongation matrices
			std::vector<SymmetricCSlCMatrix<double>> matrices; // assembled system matrices
			Index numbOfMeshLevels() { return matrices.size(); }

			boost::tuple<
				SymmetricCSlCMatrix<double>, // assembled system matrix and
				std::vector<double> // rhs vector
			> linearLagrangeSetter(
				DiffusionReactionEqn2D const & PDE,
				DirichletScalarCondition2D const & DirichletBC,
				Triangulation& Omega, // initial mesh
				Index numbOfMeshLevels
			) {
				logger.beg("mesh level 0");
					auto numbOfNodesCoarse = Omega.numbOfNodes(),
						 numbOfElementsCoarse = Omega.numbOfElements();
					auto system = linearLagrangeAssembler(PDE, Omega, DirichletBC);
					matrices.emplace_back(boost::get<0>(system));
				logger.end();
				for (Index currentMeshLevel = 1; currentMeshLevel <= numbOfMeshLevels; ++currentMeshLevel) {
					logger.beg("mesh level " + std::to_string(currentMeshLevel));
						logger.beg("refine mesh");
							Omega.refine();
						logger.end();
						logger.beg("assemble system");
							system = linearLagrangeAssembler(PDE, Omega, DirichletBC);
							matrices.emplace_back(boost::get<0>(system));
						logger.end();
						logger.beg("build pattern and fill in values of prolongation matrix");
							std::vector<Index> rowptr(Omega.numbOfNodes() + 1),
											   colind(2 * Omega.numbOfNodes());
							std::vector<double> values(colind.size(), .5); // .5 for new DOFs
							// compute row pointers, column indicies, and values for coarse DOFs
							Index i;
							for (i = 0; i < numbOfNodesCoarse; ++i) {
								rowptr[i] = colind[i] = i;
								values[i] = 1.;
							}
							rowptr[i] = i++;
							// compute row pointers for new DOFs
							for (; i < rowptr.size(); ++i) rowptr[i] = rowptr[i - 1] + 2;
							// compute column indicies for new DOFs and fix values for Dirichlet nodes
							for (Index t = 0; t < numbOfElementsCoarse; ++t) {
								auto nodesIndicies = Omega.getNodesIndicies(t);
								decltype(nodesIndicies) nodesIndiciesCoarse;
								std::array<SignedIndex, 3> neighborsIndiciesCoarse;
								for (LocalIndex j : { 0, 1, 2 }) {
									nodesIndiciesCoarse[j] = Omega.getNodesIndicies(Omega.getNeighborsIndicies(t)[j])[0];
									neighborsIndiciesCoarse[nextIndex(j)] = Omega.getNeighborsIndicies(Omega.getNeighborsIndicies(t)[j])[1];
								}
								// coarse DOFs
								for (LocalIndex j : { 0, 1, 2}) {
									auto k = nodesIndicies[j];
									//if (neighborsIndiciesCoarse[j] < 0) 
									//	values[rowptr[k]] = values[rowptr[k] + 1] = 0.;
									//else {
										colind[rowptr[k]] = nodesIndiciesCoarse[nextIndex(j)];
										colind[rowptr[k] + 1] = nodesIndiciesCoarse[nextIndex(nextIndex(j))];
									//}
								}
							}
							prolongations.emplace_back(rowptr, colind, values, numbOfNodesCoarse);
						logger.end();
						// save before the next refinement
						numbOfNodesCoarse = Omega.numbOfNodes();
						numbOfElementsCoarse = Omega.numbOfElements();
					logger.end();
				}
				return system;
			}

			std::vector<double> iteration( // for solving A.z = f
				Index meshLevel, // current mesh level
				std::vector<double> const & z_0, // initial approximation to soln
				std::vector<double> const & f // rhs
			) {
				auto numbOfCoarseDOFs = [&]() {
					return prolongations[meshLevel - 1].numbOfRows();
				};
				auto prolongate = [&](std::vector<double> const & u) {
					return prolongations[meshLevel - 1].t() * u;
				};
				auto restrict = [&](std::vector<double> const & u) {
					//return .25 * (prolongations[meshLevel - 1] * u);
					return std::vector<double>(u.begin(), u.begin() + numbOfCoarseDOFs());
				};
				if (meshLevel) { // we are not at the coarsest grid
					logger.beg("mesh level " + std::to_string(meshLevel));
						logger.beg("smooth the residual (pre-smoothing)");
							auto z = smoother(matrices[meshLevel], f, z_0);
						logger.end();
						logger.beg("restrict the residual to mesh level " + std::to_string(meshLevel - 1));
							auto d = restrict(matrices[meshLevel] * z - f);
						logger.end();
						logger.log("go to coarser grid");
						std::vector<double> e(numbOfCoarseDOFs(), 0.); // initial guess for the error on the coarser mesh level
						for (Index i = 0; i < gamma; ++i)
							e = iteration(meshLevel - 1, e, d);
						logger.beg("correct from the coarse grid");
							z -= prolongate(e);
						logger.end();
						logger.beg("smooth the residual (post-smoothing)");
							z = smoother(matrices[meshLevel], f, z);
						logger.end();
					logger.end();
					return z;
				}
				else { // we are at the coarsest grid
					logger.beg("compute exact soln");
						DenseSquareMatrix<double> B(matrices.front());
						auto z = B.GaussElimination(f);
					logger.end();
					return z;
				}
			}

		}

	}

	namespace Mixed {

		//boost::tuple<
		//	SymmetricCSlCMatrix<double>, // assembled system matrix and
		//	std::vector<double> // rhs vector
		//> 
		void
		P2P1Assembler(
			OseenProblem2D const & PDE, // (1) PDE,
			Triangulation const & Omega, // (2) discretized domain (mesh),
			DirichletVectorCondition2D const & DirichletBC // (3) and BCs that connects (1) and (2)
		) {
			// quadrature rule for master triangle
			// exact up to 4th degree polynomials
			boost::tuple<
				std::vector<Node2D>,
				std::vector<double>
			> quadratureRule( {	// nodes		
					{ .44594849091597, .44594849091597 },
					{ .44594849091597, .10810301816807 },
					{ .10810301816807, .44594849091597 },
					{ .09157621350977, .09157621350977 },
					{ .09157621350977, .81684757298046 },
					{ .81684757298046, .09157621350977 }
				}, { // weights
					0.111690794839005, 
					0.111690794839005, 
					0.111690794839005,
					0.05497587182766, 
					0.05497587182766, 
					0.05497587182766 
				}
			);
			auto& qNodes   = boost::get<0>(quadratureRule);
			auto& qWeights = boost::get<1>(quadratureRule);
			// ΔP1L shape funcs
			std::vector<ScalarField2D> pressureShapes = {
				[](Node2D const & p) { return 1. - p[0] - p[1]; },
				[](Node2D const & p) { return p[0]; },
				[](Node2D const & p) { return p[1]; },
			};
			// ΔP2L shape funcs
			std::vector<ScalarField2D> velocityShapes = {
				[](Node2D const & p) { return (-1. + p[0] + p[1]) * (-1. + 2. * p[0] + 2. * p[1]); },
				[](Node2D const & p) { return p[0] * (-1. + 2. * p[0]); },
				[](Node2D const & p) { return -1. * (1. - 2. * p[1]) * p[1]; },
				[](Node2D const & p) { return 4. * p[0] * p[1]; },
				[](Node2D const & p) { return -4. * p[1] * (-1. + p[0] + p[1]); },
				[](Node2D const & p) { return 4. * p[0] * (1. - p[0] - p[1]); }
			};
			// gradients of ″
			std::vector<VectorField2D> gradsOfVelocityShapes = {
				[](Node2D const & p) -> Node2D { return { 
						-1. + 2. * p[0] + 2. * p[1] + 2. * (-1. + p[0] + p[1]), 
						-1. + 2. * p[0] + 2. * p[1] + 2. * (-1. + p[0] + p[1]) 
				}; },
				[](Node2D const & p) -> Node2D { return { -1. + 4. * p[0], 0. }; },
				[](Node2D const & p) -> Node2D { return { 0., -1. + 4. * p[1] }; },
				[](Node2D const & p) -> Node2D { return { 4. * p[1], 4. * p[0] }; },
				[](Node2D const & p) -> Node2D { return { 
					-4. * p[1], 
					-4. * p[1] - 4. * (-1. + p[0] + p[1]) 
				}; },
				[](Node2D const & p) -> Node2D { return {
					-4. * p[0] + 4. * (1. - p[0] - p[1]), 
					-4. * p[0] 
				}; }
			};
			// values of shapes at quadrature nodes,
			// shapeVal[i][j] := value of velocityShapes[i] at qNodes[j]
			std::vector<std::vector<double>> shapeVal;
			shapeVal.reserve(velocityShapes.size());
			std::for_each(velocityShapes.begin(), velocityShapes.end(), [&](auto const & s) {
				std::vector<double> values(qNodes.size());
				std::transform(qNodes.begin(), qNodes.end(), values.begin(), [&](auto const & node) { return s(node); });
				shapeVal.emplace_back(values);
			});
			// values of gradients of shapes at quadrature nodes,
			// gradsVal[i][j] := value of ∇velocityShapes[i] at qNodes[j]
			std::vector<std::vector<Node2D>> gradsVal;
			gradsVal.reserve(velocityShapes.size());
			std::for_each(gradsOfVelocityShapes.begin(), gradsOfVelocityShapes.end(), [&](auto const & gradOfS) {
				std::vector<Node2D> values(qNodes.size());
				std::transform(qNodes.begin(), qNodes.end(), values.begin(), [&](auto const & node) { return gradOfS(node); });
				gradsVal.emplace_back(values);
			});
			// compute local mass matrix
			std::vector<double> values(qNodes.size());
			SymmetricMatrix<double> localMassMatrix(velocityShapes.size());
			for (LocalIndex i = 0; i < velocityShapes.size(); ++i)
				for (LocalIndex j = i; j < velocityShapes.size(); ++j) {
					std::transform(shapeVal[j].begin(), shapeVal[j].end(), shapeVal[i].begin(), values.begin(), std::multiplies<double>());
					localMassMatrix(i, j) = qWeights * values;
				}
			SymmetricMatrix<double> localStiffnessMatrix(velocityShapes.size());
			for (Index t = 0; t < Omega.numbOfElements(); ++t) {
				auto nodes = Omega.getElement(t);
				// jacobian–related
				auto detJ = 2. * area(Omega.getElement(t));
				DenseSquareMatrix<double> JInverseTranspose {
					{ nodes[2][1] - nodes[0][1], nodes[0][1] - nodes[1][1] },
					{ nodes[0][0] - nodes[2][0], nodes[1][0] - nodes[0][0] }
				};
				JInverseTranspose /= detJ;
				// compute local stiffness matrix
				for (LocalIndex i = 0; i < velocityShapes.size(); ++i)
					for (LocalIndex j = i; j < velocityShapes.size(); ++j) {
						std::transform(gradsVal[j].begin(), gradsVal[j].end(), gradsVal[i].begin(), values.begin(), [&](auto const & v1, auto const & v2) { return (JInverseTranspose * v1) * (JInverseTranspose * v2); });
						localStiffnessMatrix(i, j) = detJ * qWeights * values;
					}

				if (t == 8) {
					JInverseTranspose.export(logger.buf);
					logger.log();
					localStiffnessMatrix.export(logger.buf);
					logger.log();
				}

			}
		}

	}

}

//SymmetricMatrix<double> FEM::computeLocalMassMatrix(Function reactionTerm, 
//                                                       array<Node, 3>& nodes, 
//                                                       double area) {
//	SymmetricMatrix<double> m(3);
//	m(0, 0) = area * (6. * reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
//	m(0, 1) = area * (2. * reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + reactionTerm(nodes[2])) / 60.;
//	m(0, 2) = area * (2. * reactionTerm(nodes[0]) + reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
//	m(1, 1) = area * (2. * reactionTerm(nodes[0]) + 6. * reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
//	m(1, 2) = area * (reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + 2. * reactionTerm(nodes[2])) / 60.;
//	m(2, 2) = area * (2. * reactionTerm(nodes[0]) + 2. * reactionTerm(nodes[1]) + 6. * reactionTerm(nodes[2])) / 60.;
//	return m;
//}
//
//SymmetricMatrix<double> FEM::computeLocalStiffnessMatrix(Function diffusionTerm,
//                                                            array<Node, 3>& nodes,
//                                                            array<Node, 3>& middleNodes,
//                                                            double area) {
//	SymmetricMatrix<double> s(3);
//	s(0, 0) = (nodes[1].x() - nodes[2].x()) * (nodes[1].x() - nodes[2].x()) +
//	          (nodes[1].y() - nodes[2].y()) * (nodes[1].y() - nodes[2].y());
//	s(0, 1) = (nodes[0].x() - nodes[2].x()) * (nodes[2].x() - nodes[1].x()) +
//	          (nodes[0].y() - nodes[2].y()) * (nodes[2].y() - nodes[1].y());
//	s(0, 2) = (nodes[0].x() - nodes[1].x()) * (nodes[1].x() - nodes[2].x()) +
//	          (nodes[0].y() - nodes[1].y()) * (nodes[1].y() - nodes[2].y());
//	s(1, 1) = (nodes[0].x() - nodes[2].x()) * (nodes[0].x() - nodes[2].x()) +
//	          (nodes[0].y() - nodes[2].y()) * (nodes[0].y() - nodes[2].y());
//	s(1, 2) = (nodes[1].x() - nodes[0].x()) * (nodes[0].x() - nodes[2].x()) +
//	          (nodes[1].y() - nodes[0].y()) * (nodes[0].y() - nodes[2].y());
//	s(2, 2) = (nodes[0].x() - nodes[1].x()) * (nodes[0].x() - nodes[1].x()) +
//	          (nodes[0].y() - nodes[1].y()) * (nodes[0].y() - nodes[1].y());
//	for (LocalIndex i = 0; i < 3; ++i)
//		for (LocalIndex j = i; j < 3; ++j)
//			// quadratures calculated assuming diffusionTerm(x, y) 
//			// ( where (x, y) is a point from triangle spanned by nodes array )  
//			// lives in P_1, 
//			// i.e. diffusionTerm(x, y) is linear combination of {x, y, 1}:
//			/*
//			s(i, j) *= (diffusionTerm(nodes[0]) + diffusionTerm(nodes[1]) + diffusionTerm(nodes[2])) / area / 12.;
//			*/
//			// but we can do better:
//			// quadratures calculated assuming diffusionTerm(x, y) lives in P_2(ith triangle), 
//			// i.e. diffusionTerm(x, y) is linear combination of {x^2, y^2, xy, x, y, 1}:
//			s(i, j) *= (diffusionTerm(middleNodes[0]) + diffusionTerm(middleNodes[1]) + diffusionTerm(middleNodes[2])) / area / 12.;
//	return s;
//}
//
//array<double, 3> FEM::computeLocalLoadVector(Function forceTerm,
//                                             array<Node, 3>& nodes,
//                                             array<Node, 3>& middleNodes,
//                                             double area) {
//	array<double, 3> f;
//	// …assuming forceTerm(x, y) lives in P_1:
//	/*
//	return (f = {
//	2. * forceTerm(nodes[0]) + forceTerm(nodes[1]) + forceTerm(nodes[2]),
//	     forceTerm(nodes[0]) + 2. * forceTerm(nodes[1]) + forceTerm(nodes[2]),
//	     forceTerm(nodes[0]) + forceTerm(nodes[1]) + 2. * forceTerm(nodes[2])
//	}) *= area / 12.; 
//	*/
//	// …assuming forceTerm(x, y) lives in P_2(ith triangle):
//	return (f = {
//		2. * forceTerm(nodes[0]) - forceTerm(nodes[1]) - forceTerm(nodes[2]) +
//		4. * forceTerm(middleNodes[0]) + 8. * (forceTerm(middleNodes[1]) + forceTerm(middleNodes[2])),
//		2. * forceTerm(nodes[1]) - forceTerm(nodes[0]) - forceTerm(nodes[2]) +
//		4. * (2. * (forceTerm(middleNodes[0]) + forceTerm(middleNodes[2])) + forceTerm(middleNodes[1])),
//		2. * (forceTerm(nodes[2]) + 4. * (forceTerm(middleNodes[0]) + forceTerm(middleNodes[1])) + 2. * forceTerm(middleNodes[2])) -
//		forceTerm(nodes[0]) - forceTerm(nodes[1])
//	}) *= area / 60.;
//}
//
//
//SymmetricMatrix<double> FEM::computeLocalRobinMatrix(BoundaryConditions const & BCs,
//                                                        array<Node, 2>& nodes,
//                                                        double length) {
//	SymmetricMatrix<double> r(2);
//	r(0, 0) = length * ( 3. * BCs.RobinCoefficient(nodes[0]) +      BCs.RobinCoefficient(nodes[1]) ) / 12.;
//	r(0, 1) = length * (      BCs.RobinCoefficient(nodes[0]) +      BCs.RobinCoefficient(nodes[1]) ) / 12.;
//	r(1, 1) = length * (      BCs.RobinCoefficient(nodes[0]) + 3. * BCs.RobinCoefficient(nodes[1]) ) / 12.;
//	return r;
//}
//
//array<double, 2> FEM::computeLocalRobinVector(BoundaryConditions const & BCs,
//                                              array<Node, 2>& nodes,
//                                              double length) {
//	array<double, 2> r;
//	return (r = {
//		4. * BCs.NeumannValue(nodes[0]) + 2. * BCs.NeumannValue(nodes[1]) +
//		BCs.DirichletCondition(nodes[0]) * ( 3. * BCs.RobinCoefficient(nodes[0]) + BCs.RobinCoefficient(nodes[1]) ) +
//		BCs.DirichletCondition(nodes[1]) * ( BCs.RobinCoefficient(nodes[0]) + BCs.RobinCoefficient(nodes[1]) ),
//		2. * BCs.NeumannValue(nodes[0]) + 4. * BCs.NeumannValue(nodes[1]) +
//		BCs.DirichletCondition(nodes[0]) * ( BCs.RobinCoefficient(nodes[0]) + BCs.RobinCoefficient(nodes[1]) ) +
//		BCs.DirichletCondition(nodes[1]) * ( BCs.RobinCoefficient(nodes[0]) + 3. * BCs.RobinCoefficient(nodes[1]) )
//	}) *= length / 12.;
//}
//
//vector<double> FEM::constructVector(Function u, Triangulation& Omega) {
//	vector<double> uVec(Omega.numbOfNodes());
//	for (size_t i = 0; i < uVec.size(); ++i)
//		uVec[i] = u(Omega.getNode(i));
//	return uVec;
//}
//
//vector<size_t> FEM::computeBoundaryNodes(Triangulation& Omega, Boundary& bndry, Predicate validAt) {
//	vector<size_t> res;
//	for (auto const & edge : bndry) 
//		if (validAt(Omega.getNode(edge[0]))) res.push_back(edge[0]);
//	return res;
//}