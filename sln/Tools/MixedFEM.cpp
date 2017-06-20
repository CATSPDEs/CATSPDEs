#include "MixedFEM.hpp"
// logger
#include "SingletonLogger.hpp"
// quadrature rule 
#include "MasterTriangleQuadratureRule.hpp"
#include "MasterSegmentQuadratureRule.hpp"
// linear shapes for transformation from master to physical triangle
#include "Triangle_P1_Lagrange.hpp"
// for local matrices
#include "SymmetricMatrix.hpp"
#include "DenseMatrix.hpp"

namespace FEM {

	namespace Mixed {

		boost::tuple<
			CSlCMatrix<double>, // diffusion + convection + reaction matrix
			CSCMatrix<double>, CSCMatrix<double>, // divergence matrices 
			std::vector<double> // rhs vector
		>
		assembleSystem(
			OseenProblem2D const & PDE, // (1) PDE,
			Triangulation const & Omega, // (2) discretized domain (mesh), and BCs that connects (1) and (2):
			VectorBoundaryCondition2D const & NeumannBC, // natural (Neumann) BC,
			VectorBoundaryCondition2D const & DirichletBC, // strong (Dirichlet) BC
			TriangularScalarFiniteElement const & velocityFE, // for each velocity component
			TriangularScalarFiniteElement const & pressureFE, // for pressure
			boost::optional<Index&> activeElementIndex
		) {
			// logger
			auto& logger = SingletonLogger::instance();
			// master element
			Triangle2D master {
				{ { 0., 0. },{ 1., 0. },{ 0., 1. } }
			};
			// r[i] : [-1, 1] —> ith edge of master triangle,
			// T(r[i]) (see below) := mapping from ith edge of master triangle to ith edge of phisical ″
			std::vector<Curve2D> r {
				[](double const & s) -> Node2D { return { .5 * (1. - s), .5 * (1. + s) }; },
				[](double const & s) -> Node2D { return { 0., .5 * (1. - s) }; },
				[](double const & s) -> Node2D { return { .5 * (1. + s), 0. }; }
			};
			// master shapes for mapping from master element to physical ″
			auto lagrangeMasterShapes = Triangle_P1_Lagrange::instance().getShapesOf(master);
			// shapes for velocity components and pressure
			auto pressureMasterShapes = pressureFE.getShapesOf(master);
			auto velocityMasterShapes = velocityFE.getShapesOf(master);
			auto velocityMasterSGrads = velocityFE.getSGradsOf(master);
			// Convection matrix integrand has max polynomial degree
			// so we need it in order to choose quadrature rule
			LocalIndex deg = ceil(velocityFE.deg() + (velocityFE.deg() - 1.) + velocityFE.deg());
			logger.buf << "polynomial degree for Gaussian quadrature: " << deg;
			logger.log();
			// set up quadrature rule
			logger.beg("precompute values of shapes at quadrature nodes");
				// quadrature rule for master triangle
				auto& qRuleTriangle = MasterTriangleQuadratureRule::instance();
				auto qNodesTriangle = qRuleTriangle.getQuadratureNodes(deg);
				// ″ for master segment
				auto& qRuleSegment = MasterSegmentQuadratureRule::instance();
				auto qNodesSegment = qRuleSegment.getQuadratureNodes(deg);
				// we will also need images of velocity shapes on edges (for quadratures due to Neumann BCs)
				// so we should compute quadrature nodes on edges, too
				decltype(qNodesTriangle) qBoundaryNodesTriangle;
				qBoundaryNodesTriangle.reserve(3 * qNodesSegment.size());
				for (LocalIndex i : { 0, 1, 2 }) // for each edge of master triangle
					for (auto const & p : qNodesSegment) // and quadrature node (1D) on master segment,
						qBoundaryNodesTriangle.emplace_back(r[i](p)); // compute corresponding quadrature node on master triangle (2D) using curves r
				// so now we are ready to precompute images of q nodes
				for (auto& s : lagrangeMasterShapes) s.saveImagesOf(qNodesTriangle).saveImagesOf(qBoundaryNodesTriangle);
				for (auto& s : velocityMasterShapes) s.saveImagesOf(qNodesTriangle).saveImagesOf(qBoundaryNodesTriangle);
				for (auto& g : velocityMasterSGrads) g.saveImagesOf(qNodesTriangle);
				for (auto& s : pressureMasterShapes) s.saveImagesOf(qNodesTriangle);
			logger.end();
			// local matrices and vectors
			SymmetricMatrix<double> localMassMatrixSkeleton(velocityMasterShapes.size()),
			                        localMassMatrix(velocityMasterShapes.size()),
			                        localStiffnessMatrix(velocityMasterShapes.size());
			DenseMatrix<double>     localConvectionMatrix(velocityMasterShapes.size()),
			                        localDivergenceMatrix1(pressureMasterShapes.size(), velocityMasterShapes.size()),
			                        localDivergenceMatrix2(pressureMasterShapes.size(), velocityMasterShapes.size());
			std::vector<double>     localLoadVector1(velocityMasterShapes.size()),
			                        localLoadVector2(velocityMasterShapes.size()),
			                        localContVector(pressureMasterShapes.size());
			// compute skeleton of local mass matrix
			for (LocalIndex i = 0; i < velocityMasterShapes.size(); ++i)
				for (LocalIndex j = i; j < velocityMasterShapes.size(); ++j)
					localMassMatrixSkeleton(i, j) = qRuleTriangle.computeQuadrature([&](Node2D const & p) {
						return velocityMasterShapes[j](p) * velocityMasterShapes[i](p);
					}, deg);
			// for Dirichlet BCs
			std::array<Index2Value<double>, 2> ind2val;
			std::array<std::vector<LocalIndex>, 3> velocityBndryDOFsLocalIndicies;
			for (LocalIndex i : {0, 1, 2}) velocityBndryDOFsLocalIndicies[i] = velocityFE.getBndryDOFsLocalIndicies(Omega, i);
			// current element and its middle nodes
			Triangle2D enodes, mnodes;
			// ribs of the current element
			std::array<Segment2D, 3> ribs;
			std::vector<Index>  velocityDOFsNumn;
			std::vector<Node2D> velocityDOFsNodes;
			// mapping from master element to physical ″
			VectorField2D T = [&](Node2D const & masterNode) {
				Node2D physicalNode{ { 0., 0. } };
				for (LocalIndex i = 0; i < 3; ++i) {
					physicalNode[0] += enodes[i][0] * lagrangeMasterShapes[i](masterNode);
					physicalNode[1] += enodes[i][1] * lagrangeMasterShapes[i](masterNode);
				}
				return physicalNode;
			};
			// jacobian–related, J := jacobian matrix of the mapping T
			double detJ;
			DenseMatrix<double> JInverseTranspose;
			// final system
			auto n = velocityFE.numbOfDOFs(Omega);
			auto m = pressureFE.numbOfDOFs(Omega);
			
			CSlCMatrix<double> A11(n);
			//CSCMatrix<double> A11(n, n);
			
			CSCMatrix<double> B1(m, n), B2(m, n);
			std::vector<double> f(2 * n + m);
			// build pattern
			logger.beg("build pattern of matrices");
				logger.beg("build pattern of A11 block");
					A11.generatePatternFrom(
						createDOFsConnectivityList(Omega, velocityFE)
					);
					//A11.generatePatternFrom(
					//	createDOFsConnectivityList(Omega, velocityFE, velocityFE, [](Index, Index) { return true; })
					//);
					logger.buf << "A11 size    = " << n << '\n'
					           << "numb of nnz = " << A11.nnz();
					logger.log();
				logger.end();
				logger.beg("build pattern of B1 and B2 blocks");
					B1.generatePatternFrom(createDOFsConnectivityList(
						Omega, velocityFE, pressureFE
					));
					B2 = B1;
					logger.buf << "B1/B2 sizes = " << m << " * " << n << '\n'
					           << "numb of nnz = " << B1.nnz();
					logger.log();
				logger.end();
				logger.buf << "overall size = " << 2 * n + m;
				logger.log();
			logger.end();
			// assemble
			logger.beg("assemble system matrix and rhs vector");
				// index of the active element
				Index activeElementIndexScoped;
				Index& t = activeElementIndex.value_or(activeElementIndexScoped);
				for (t = 0; t < Omega.numbOfElements(); ++t) {
					enodes = Omega.getElement(t);
					mnodes = midNodes(enodes);
					ribs = ribsOf(enodes);
					velocityDOFsNumn = velocityFE.getDOFsNumeration(Omega, t);
					velocityDOFsNodes = velocityFE.getDOFsNodes(Omega, t);
					detJ = 2. * area(Omega.getElement(t));
					JInverseTranspose = {
						{ enodes[2][1] - enodes[0][1], enodes[0][1] - enodes[1][1] },
						{ enodes[0][0] - enodes[2][0], enodes[1][0] - enodes[0][0] }
					};
					JInverseTranspose /= detJ;
					// compute local mass matrix
					localMassMatrix = PDE.massCoef() * detJ * localMassMatrixSkeleton;
					// compute local stiffness matrix
					for (LocalIndex i = 0; i < velocityMasterShapes.size(); ++i)
						for (LocalIndex j = i; j < velocityMasterShapes.size(); ++j)
							localStiffnessMatrix(i, j) = PDE.inverseReynoldsNumber() * detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return (JInverseTranspose * velocityMasterSGrads[j](p)) * (JInverseTranspose * velocityMasterSGrads[i](p));
							}, deg);
					// compute local convection matrix
					for (LocalIndex i = 0; i < velocityMasterShapes.size(); ++i)
						for (LocalIndex j = 0; j < velocityMasterShapes.size(); ++j) 
							localConvectionMatrix(i, j) = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return (PDE.windField(T(p)) * (JInverseTranspose * velocityMasterSGrads[j](p))) * velocityMasterShapes[i](p);
							}, deg);
					// compute local divergence matrix
					for (LocalIndex i = 0; i < pressureMasterShapes.size(); ++i)
						for (LocalIndex j = 0; j < velocityMasterShapes.size(); ++j) {
							localDivergenceMatrix1(i, j) = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return (JInverseTranspose.getRow(0) * velocityMasterSGrads[j](p)) * pressureMasterShapes[i](p);
							}, deg);
							localDivergenceMatrix2(i, j) = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return (JInverseTranspose.getRow(1) * velocityMasterSGrads[j](p)) * pressureMasterShapes[i](p);
							}, deg);
						}
					// compute local load vector
					for (LocalIndex i = 0; i < velocityMasterShapes.size(); ++i) {
						localLoadVector1[i] = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
							return PDE.forceTerm(T(p))[0] * velocityMasterShapes[i](p);
						}, deg);
						localLoadVector2[i] = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
							return PDE.forceTerm(T(p))[1] * velocityMasterShapes[i](p);
						}, deg);
					}
					// compute local “continuity” vector
					for (LocalIndex i = 0; i < pressureMasterShapes.size(); ++i)
						localContVector[i] = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
							return PDE.continuityTerm(T(p)) * pressureMasterShapes[i](p);
						}, deg);
					// compute quadratures over boundary
					for (LocalIndex i : { 0, 1, 2 }) // loop over ribs
						if (Omega.getNeighborsIndicies(t)[i] >= 0) continue;
						else if (NeumannBC.shouldBeEnforcedAt(mnodes[i]))
							for (LocalIndex j : velocityBndryDOFsLocalIndicies[i]) {
								localLoadVector1[j] += .5 * length(ribs[i]) * qRuleSegment.computeQuadrature([&](double const & s) {
									return NeumannBC(T(r[i](s)))[0] * velocityMasterShapes[j](r[i](s));
								}, deg);
								localLoadVector2[j] += .5 * length(ribs[i]) * qRuleSegment.computeQuadrature([&](double const & s) {
									return NeumannBC(T(r[i](s)))[1] * velocityMasterShapes[j](r[i](s));
								}, deg);
							}
						else if (DirichletBC.shouldBeEnforcedAt(mnodes[i]))
							for (LocalIndex j : velocityBndryDOFsLocalIndicies[i]) {
								ind2val[0][velocityDOFsNumn[j]] = DirichletBC(velocityDOFsNodes[j])[0];
								ind2val[1][velocityDOFsNumn[j]] = DirichletBC(velocityDOFsNodes[j])[1];
							}
						else logger.wrn("no BCs were prescribed; hom. Neumann was assumed");
					// to global system
					auto l2gVel = velocityFE.getDOFsNumeration(Omega, t);
					auto l2gPre = pressureFE.getDOFsNumeration(Omega, t);
					for (LocalIndex i = 0; i < velocityMasterShapes.size(); ++i) {
						for (LocalIndex j = 0; j < velocityMasterShapes.size(); ++j)
							A11(l2gVel[i], l2gVel[j]) += localMassMatrix(i, j) + localStiffnessMatrix(i, j) + localConvectionMatrix(i, j);
						f[l2gVel[i]    ] += localLoadVector1[i];
						f[l2gVel[i] + n] += localLoadVector2[i];
					}
					for (LocalIndex i = 0; i < pressureMasterShapes.size(); ++i) {
						for (LocalIndex j = 0; j < velocityMasterShapes.size(); ++j) {
							B1(l2gPre[i], l2gVel[j]) -= localDivergenceMatrix1(i, j);
							B2(l2gPre[i], l2gVel[j]) -= localDivergenceMatrix2(i, j);
						}
						f[2 * n + l2gPre[i]] += localContVector[i];
					}
				}
			logger.end();
			if (ind2val[0].size()) {
				logger.beg("enforce no-slip BCs");
					auto A22 = A11;
					A11.enforceDirichletBCs(ind2val[0], f.data());
					A22.enforceDirichletBCs(ind2val[1], f.data() + n);
					B1.enforceDirichletBCs(ind2val[0], f.data() + 2 * n);
					B2.enforceDirichletBCs(ind2val[1], f.data() + 2 * n);
				logger.end();
			}
			return { A11, B1, B2, f };
		}

	}

}