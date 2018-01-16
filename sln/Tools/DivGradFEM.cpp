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
// in order to store computed images of quadrature nodes
#include "SmartMapping.hpp"

namespace FEM {

	namespace DivGrad {

		boost::tuple<
			CSlCMatrix<double>, // system matrix
			std::vector<double> // rhs vector
		> assembleSystem(
			ConvectionDiffusionEqn2D const & PDE,
			Triangulation const & Omega,
			ScalarBoundaryCondition2D const & RobinBC,
			ScalarBoundaryCondition2D const & DirichletBC,
			TriangularScalarFiniteElement const & FE,
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
			auto lagrangeMasterShapes = smartify<2, 1>(Triangle_P1_Lagrange::instance().getShapesOf(master));
			// shapes of FE
			auto masterShapes = smartify<2, 1>(FE.getShapesOf(master));
			auto masterSGrads = smartify<2, 2>(FE.getSGradsOf(master));
			auto l = masterShapes.size();
			// set up quadrature rules data
			LocalIndex deg = ceil(3. * FE.deg());
			// quadrature rule for master triangle
			auto& qRuleTriangle = MasterTriangleQuadratureRule::instance();
			// ″ for master segment
			auto& qRuleSegment = MasterSegmentQuadratureRule::instance();
			Index capacity;
			logger.buf
				<< "polynomial degree for Gaussian quadrature: " << deg << '\n'
				<< "numb of q nodes (per element): " << qRuleTriangle.getQuadratureNodes(deg).size() << '\n'
				<< "numb of q nodes (per face): " << qRuleSegment.getQuadratureNodes(deg).size() << '\n'
				<< "=> numb of values to be saved (per shape func): " << (capacity = qRuleTriangle.getQuadratureNodes(deg).size() + 3 * qRuleSegment.getQuadratureNodes(deg).size());
			logger.log();
			for (auto& s : lagrangeMasterShapes) s.capacity() = capacity;
			for (auto& s : masterShapes) s.capacity() = capacity;
			for (auto& g : masterSGrads) g.capacity() = capacity;
			// local matrices
			SymmetricMatrix<double> localMassMatrix(l),
			                        localStiffnessMatrix(l);
			DenseMatrix<double>     localConvectionMatrix(l);
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
			// mapping from master element to physical ″
			VectorField2D T = [&](Node2D const & masterNode) {
				Node2D physicalNode { { 0., 0. } };
				for (LocalIndex i = 0; i < 3; ++i) {
					physicalNode[0] += enodes[i][0] * lagrangeMasterShapes[i](masterNode);
					physicalNode[1] += enodes[i][1] * lagrangeMasterShapes[i](masterNode);
				}
				return physicalNode;
			};
			// jacobian–related, J := jacobian matrix of the mapping T
			double detJ;
			DenseMatrix<double> JInverseTranspose;
			// data structures for final linear system A.xi = b:
			logger.beg("generate matrix pattern");
				auto n = FE.numbOfDOFs(Omega);
				CSlCMatrix<double> A(n); // system matrix
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
					detJ = 2. * area(Omega.getElement(t));
					JInverseTranspose = {
						{ enodes[2][1] - enodes[0][1], enodes[0][1] - enodes[1][1] },
						{ enodes[0][0] - enodes[2][0], enodes[1][0] - enodes[0][0] }
					};
					JInverseTranspose /= detJ;
					// compute local stiffness matrix
					for (LocalIndex i = 0; i < l; ++i) {
						for (LocalIndex j = i; j < l; ++j) {
							localStiffnessMatrix(i, j) = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return PDE.diffusion(T(p)) * (JInverseTranspose * masterSGrads[j](p)) * (JInverseTranspose * masterSGrads[i](p));
							}, deg);
							localMassMatrix(i, j) = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return PDE.reaction(T(p)) * masterShapes[j](p) * masterShapes[i](p);
							}, deg);
						}
						// compute local convection matrix
						for (LocalIndex j = 0; j < l; ++j)
							localConvectionMatrix(i, j) = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
								return (PDE.convection(T(p)) * (JInverseTranspose * masterSGrads[j](p))) * masterShapes[i](p);
							}, deg);
						localLoadVector[i] = detJ * qRuleTriangle.computeQuadrature([&](Node2D const & p) {
							return PDE.force(T(p)) * masterShapes[i](p);
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
						for (LocalIndex j = 0; j < l; ++j)
							A(DOFsNumn[i], DOFsNumn[j]) += localMassMatrix(i, j) + localStiffnessMatrix(i, j) + localConvectionMatrix(i, j);
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