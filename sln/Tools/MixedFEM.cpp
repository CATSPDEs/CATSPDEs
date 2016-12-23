#include "MixedFEM.hpp"
// logger
#include "SingletonLogger.hpp"
// quadrature rule 
#include "MasterTriangleQuadratureRule.hpp"
// linear shapes for transformation from master to physical triangle
#include "Triangle_P1_Lagrange.hpp"
// for local matrices
#include "SymmetricMatrix.hpp"
#include "DenseMatrix.hpp"

namespace FEM {

	namespace Mixed {

		//boost::tuple<
		//	SymmetricCSlCMatrix<double>, // assembled system matrix and
		//	std::vector<double> // rhs vector
		//> 
		void
		assembleSystem(
			OseenProblem2D const & PDE, // (1) PDE,
			Triangulation const & Omega, // (2) discretized domain (mesh), and BCs that connects (1) and (2):
			VectorBoundaryCondition2D const & DirichletBC, // strong (Dirichlet) BC, 
			VectorBoundaryCondition2D const & NeumannBC, // natural (Neumann) BC
			TriangularScalarFiniteElement const & velocityFE, // for each velocity component
			TriangularScalarFiniteElement const & pressureFE// for pressure
		) {
			// logger
			auto& logger = SingletonLogger::instance();
			// master element
			Triangle master {
				{ { 0., 0. }, { 1., 0. }, { 0., 1. } }
			};
			// master shapes for mapping from master element to physical ″
			auto lagrangeShapes = Triangle_P1_Lagrange::instance().getShapesOf(master);
			// shapes for velocity components and pressure
			auto pressureShapes = pressureFE.getShapesOf(master);
			auto velocityShapes = velocityFE.getShapesOf(master);
			auto velocitySGrads = velocityFE.getSGradsOf(master);
			// Newton matrix integrand has max polynomial degree
			// so we need it in order to choose quadrature rule
			LocalIndex deg = ceil(velocityFE.deg() + (velocityFE.deg() - 1.) + velocityFE.deg());
			logger.buf << "polynomial degree for Gaussian quadrature: " << deg;
			logger.log();
			// set up quadrature rule
			logger.beg("precompute values of shapes at quadrature nodes");
				auto& qRule = MasterTriangleQuadratureRule::instance();
				auto qNodes = qRule.getQuadratureNodes(deg);
				for (auto& s : lagrangeShapes) s.saveImagesOf(qNodes);
				for (auto& s : velocityShapes) s.saveImagesOf(qNodes);
				for (auto& g : velocitySGrads) g.saveImagesOf(qNodes);
				for (auto& s : pressureShapes) s.saveImagesOf(qNodes);
			logger.end();
			// local matrices
			SymmetricMatrix<double> localMassMatrixSkeleton(velocityShapes.size()),
			                        localMassMatrix(velocityShapes.size()),
			                        localStiffnessMatrix(velocityShapes.size());
			DenseMatrix<double>     localNewtonMatrix(velocityShapes.size());
			// compute skeleton of local mass matrix
			for (LocalIndex i = 0; i < velocityShapes.size(); ++i)
				for (LocalIndex j = i; j < velocityShapes.size(); ++j)
					localMassMatrixSkeleton(i, j) = qRule.computeQuadrature([&](Node2D const & p) {
						return velocityShapes[j](p) * velocityShapes[i](p);
					}, deg);
			// current element
			Triangle element; 
			// mapping from master element to physical ″
			VectorField2D T = [&](Node2D const & masterNode) { 
				Node2D physicalNode { { 0., 0. } };
				for (LocalIndex i = 0; i < 3; ++i) {
					physicalNode[0] += element[i][0] * lagrangeShapes[i](masterNode);
					physicalNode[1] += element[i][1] * lagrangeShapes[i](masterNode);
				}
				return physicalNode;
			}; 
			// assemble
			logger.beg("assemble system matrix and rhs vector");
				for (Index t = 0; t < Omega.numbOfElements(); ++t) {
					element = Omega.getElement(t);
					// jacobian–related, J := jacobian of the mapping T
					auto detJ = 2. * area(Omega.getElement(t));
					DenseMatrix<double> JInverseTranspose {
						{ element[2][1] - element[0][1], element[0][1] - element[1][1] },
						{ element[0][0] - element[2][0], element[1][0] - element[0][0] }
					};
					JInverseTranspose /= detJ;
					// compute local mass matrix
					localMassMatrix = PDE.massTerm() * detJ * localMassMatrixSkeleton;
					// compute local stiffness matrix
					for (LocalIndex i = 0; i < velocityShapes.size(); ++i)
						for (LocalIndex j = i; j < velocityShapes.size(); ++j)
							localStiffnessMatrix(i, j) = PDE.inverseReynoldsNumber() * detJ * qRule.computeQuadrature([&](Node2D const & p) {
								return (JInverseTranspose * velocitySGrads[j](p)) * (JInverseTranspose * velocitySGrads[i](p));
							}, deg);
					// compute local Newton matrix
					for (LocalIndex i = 0; i < velocityShapes.size(); ++i)
						for (LocalIndex j = 0; j < velocityShapes.size(); ++j) 
							localNewtonMatrix(i, j) = detJ * qRule.computeQuadrature([&](Node2D const & p) {
								return (PDE.NewtonTerm(T(p)) * (JInverseTranspose * velocitySGrads[j](p))) * velocityShapes[i](p);
							}, deg);
					// temp
					if (t == 12) {
						localMassMatrix.export(logger.buf);
						localStiffnessMatrix.export(logger.buf);
						localNewtonMatrix.export(logger.buf);
						logger.log();
					}
				}
			logger.end();
		}

	}

}