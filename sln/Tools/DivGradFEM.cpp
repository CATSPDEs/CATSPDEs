#include "DivGradFEM.hpp"
// logger
#include "SingletonLogger.hpp"
// for local matrices
#include "SymmetricMatrix.hpp"
#include "DenseMatrix.hpp"
// quadrature rule 
#include "MasterTriangleQuadratureRule.hpp"
#include "MasterSegmentQuadratureRule.hpp"
// linear shapes for transformation from master to physical triangle
#include "Triangle_P1_Lagrange.hpp"

namespace FEM {

	namespace DivGrad {

		boost::tuple<
			SymmetricCSlCMatrix<double>, // system matrix
			std::vector<double> // rhs vector
		> assembleSystem(
			DiffusionReactionEqn2D const & PDE,
			Triangulation const & Omega,
			ScalarBoundaryCondition2D const & RobinBC,
			ScalarBoundaryCondition2D const & DirichletBC,
			TriangularScalarFiniteElement const & FE,
			boost::optional<TriangularScalarFiniteElement const &> T_FE,
			boost::optional<Index&> activeElementIndex
		) {
			// logger
			auto& logger = SingletonLogger::instance();
			// master element
			Triangle2D master {
				{ { 0., 0. }, { 1., 0. }, { 0., 1. } }
			};
			// r[i] : [-1, 1] —> ith edge of master triangle,
			// T(r[i]) (see below) := mapping from ith edge of master triangle to ith edge of phisical ″
			std::vector<Curve2D> r {
				[](double const & s) -> Node2D { return { .5 * (1. - s), .5 * (1. + s) }; },
				[](double const & s) -> Node2D { return { 0., .5 * (1. - s) }; },
				[](double const & s) -> Node2D { return { .5 * (1. + s), 0. }; }
			};
			// shapes for transformation T
			auto T_masterShapes = T_FE.value_or(FE).getShapesOf(master);
			auto T_masterSGrads = T_FE.value_or(FE).getSGradsOf(master);
			auto T_l = T_masterShapes.size();
				// T_FE.value_or(Triangle_P1_Lagrange::instance()).getShapesOf(master);
			// shapes of FE
			auto masterShapes = FE.getShapesOf(master);
			auto masterSGrads = FE.getSGradsOf(master);
			auto l = masterShapes.size();
			LocalIndex deg = ceil(3. * FE.deg());
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
				// we will also need images of shapes on edges (for quadratures due to Neumann BCs)
				// so we should compute quadrature nodes on edges, too
				decltype(qNodesTriangle) qBoundaryNodesTriangle;
				qBoundaryNodesTriangle.reserve(3 * qNodesSegment.size());
				for (LocalIndex i : { 0, 1, 2 }) // for each edge of master triangle
					for (auto const & p : qNodesSegment) // and quadrature node (1D) on master segment,
						qBoundaryNodesTriangle.emplace_back(r[i](p)); // compute corresponding quadrature node on master triangle (2D) using curves r
				// so now we are ready to precompute images of q nodes
				for (auto& s : T_masterShapes) s.saveImagesOf(qNodesTriangle).saveImagesOf(qBoundaryNodesTriangle);
				for (auto& s : masterShapes) s.saveImagesOf(qNodesTriangle).saveImagesOf(qBoundaryNodesTriangle);
				for (auto& g : masterSGrads) g.saveImagesOf(qNodesTriangle);
			logger.end();
			// local matrices
			SymmetricMatrix<double> localMassMatrix(l),
			                        localStiffnessMatrix(l);
			std::vector<double>     localLoadVector(l);
			// for Dirichlet BCs
			Index2Value<double> ind2val;
			std::array<std::vector<LocalIndex>, 3> bndryDOFsLocalIndicies;
			for (LocalIndex i : {0, 1, 2}) bndryDOFsLocalIndicies[i] = FE.getBndryDOFsLocalIndicies(Omega, i);
			// current element and its middle nodes
			Triangle2D enodes, mnodes;
			// ribs of the current element
			std::array<Segment2D, 3> ribs;
			std::vector<Index> DOFsNumn;		  
			std::vector<Node2D> DOFsNodes;	
			std::vector<Node2D> T_DOFsNodes;

			// mapping from master element to physical ″
			auto T = [&](Node2D const & p)->Node2D {
				std::vector<double> sValues(T_l);
				std::vector<double> X(T_l), Y(T_l);
				for (LocalIndex i = 0; i < T_l; i++)
				{
					X[i] = T_DOFsNodes[i][0];
					Y[i] = T_DOFsNodes[i][1];
				}
				std::transform(T_masterShapes.begin(), T_masterShapes.end(), sValues.begin(), [&](auto const & s) {
					return s(p);
				});
				return{ X*sValues ,Y*sValues };
				/*Node2D physicalNode { { 0., 0. } };
				for (LocalIndex i = 0; i < T_masterShapes.size(); ++i) {
					physicalNode[0] += enodes[i][0] * T_masterShapes[i](masterNode);
					physicalNode[1] += enodes[i][1] * T_masterShapes[i](masterNode);
				}
				return physicalNode;*/
			};
			auto gradT = [&](Node2D const & p)->DenseMatrix<double> {
				
				std::vector<Node2D> gValues(T_l);
				std::vector<double> X(T_l), Y(T_l);
				for (LocalIndex i = 0; i < T_l; i++)
				{
					X[i] = T_DOFsNodes[i][0];
					Y[i] = T_DOFsNodes[i][1];
				}

				std::transform(T_masterSGrads.begin(), T_masterSGrads.end(), gValues.begin(), [&](auto const & g) {
					return g(p);
				});

				return
				{
					{(X*gValues)[0],(X*gValues)[1] },
					{ (Y*gValues)[0],(Y*gValues)[1] }
				};


			};
			// data structures for final linear system A.xi = b:
			logger.beg("generate matrix pattern");
				auto n = FE.numbOfDOFs(Omega);
				SymmetricCSlCMatrix<double> A(n); // system matrix
				std::vector<double> b(n); // load vector
				A.generatePatternFrom(createDOFsConnectivityList(Omega, FE));
				logger.buf << "system size = " << n << '\n'
				           << "numb of nnz = " << A.nnz();
				logger.log();
			logger.end();
			logger.beg("assemble system matrix and rhs vector");
				// index of the active element
				Index activeElementIndexScoped;
				Index& t = activeElementIndex.value_or(activeElementIndexScoped);
				for (t = 0; t < Omega.numbOfElements(); ++t) {
					enodes = Omega.getElement(t);
					mnodes = midNodes(enodes);
					ribs = ribsOf(enodes);
					DOFsNumn = FE.getDOFsNumeration(Omega, t);
					DOFsNodes = FE.getDOFsNodes(Omega, t);
					T_DOFsNodes= T_FE.value_or(FE).getDOFsNodes(Omega, t);

					/*detJ = 2. * area(Omega.getElement(t));
					JInverseTranspose = {
						{ enodes[2][1] - enodes[0][1], enodes[0][1] - enodes[1][1] },
						{ enodes[0][0] - enodes[2][0], enodes[1][0] - enodes[0][0] }
					};*/
					//JInverseTranspose /= detJ;


					// compute local stiffness matrix
					for (LocalIndex i = 0; i < l; ++i) {
						for (LocalIndex j = i; j < l; ++j) {
							
							//localStiffnessMatrix(i, j) = detJ * qRuleTriangle.computeQuadrature(stiffnessMatrixIntegrand[i][j], deg);
							//localMassMatrix(i, j) = detJ * qRuleTriangle.computeQuadrature(massMatrixIntegrand[i][j], deg);
							
							localStiffnessMatrix(i, j) = qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return PDE.diffusionTerm(T(p)) * (gradT(p).inv().t() * masterSGrads[j](p)) * (gradT(p).inv().t() * masterSGrads[i](p)) * fabs(gradT(p).det());
							}, deg);
							localMassMatrix(i, j) =qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return PDE.reactionTerm(T(p)) * masterShapes[j](p) * masterShapes[i](p)*fabs(gradT(p).det());
							}, deg);
						}
						
						//localLoadVector[i] = detJ * qRuleTriangle.computeQuadrature(loadVectorIntegrand[i], deg);
						
						localLoadVector[i] = qRuleTriangle.computeQuadrature([&](Node2D const & p) {
							return PDE.forceTerm(T(p)) * masterShapes[i](p)*fabs(gradT(p).det());
						}, deg);
					}
					// assemble boundary data
					for (LocalIndex i : { 0, 1, 2 }) // loop over ribs
						if (Omega.getNeighborsIndicies(t)[i] >= 0) continue;
						else if (RobinBC.shouldBeEnforcedAt(mnodes[i])) 
							for (LocalIndex j : bndryDOFsLocalIndicies[i]) {
								for (LocalIndex k : bndryDOFsLocalIndicies[i]) 
									if (k >= j)
										localMassMatrix(j, k) += .5 * length(ribs[i]) * qRuleSegment.computeQuadrature([&](double const & s) {
											return RobinBC(T(r[i](s))) * masterShapes[k](r[i](s)) * masterShapes[j](r[i](s));
										}, deg);
								localLoadVector[j] += .5 * length(ribs[i]) * qRuleSegment.computeQuadrature([&](double const & s) {
									return RobinBC(T(r[i](s)), 1) * masterShapes[j](r[i](s));
								}, deg);
							}
						else if (DirichletBC.shouldBeEnforcedAt(mnodes[i]))
							for (LocalIndex j : bndryDOFsLocalIndicies[i])
								ind2val[DOFsNumn[j]] = DirichletBC(DOFsNodes[j]);
						else logger.wrn("no BCs were prescribed; hom. Neumann was assumed");
					// to global system
					for (LocalIndex i = 0; i < l; ++i) {
						for (LocalIndex j = i; j < l; ++j)
							A(DOFsNumn[i], DOFsNumn[j]) += localMassMatrix(i, j) + localStiffnessMatrix(i, j);
						b[DOFsNumn[i]] += localLoadVector[i];
					}
				}
			logger.end();
			if (ind2val.size()) {
				logger.beg("enforce strong BCs");
					A.enforceDirichletBCs(ind2val, b.data());
				logger.end();
			}
			return { A, b };
		}

	}

}