#include "L2Projection.hpp"
// logger
#include "SingletonLogger.hpp"
// quadrature rule
#include "MasterTriangleQuadratureRule.hpp"
// for local matrices
#include "SymmetricMatrix.hpp"
#include "DenseMatrix.hpp"
// linear shapes for transformation from master to physical triangle
#include "Triangle_P1_Lagrange.hpp"

/*
	Žilyakov Alexander, Jan 2017
*/

namespace FEM {

	boost::tuple<
		SymmetricCSlCMatrix<double>, // system mass matrix
		CSCMatrix<double> // rhs mass matrix
	>
	L2ProjectionAssembler(
		Triangulation const & fMesh, // fine and
		Triangulation const & cMesh, // coarse meshes
		TriangularScalarFiniteElement const & FE // FE to define interpolants
	) {
		// logger
		auto& logger = SingletonLogger::instance();
		// resulting matrices
		auto m = FE.numbOfDOFs(cMesh), n = FE.numbOfDOFs(fMesh);
		SymmetricCSlCMatrix<double> cM(m); // coarse mass matrix
		CSCMatrix<double> cfM(m, n); // coarse–fine mass matrix
		// so let u := dofs of fine FE–interpolant, v := ″ of coarse,
		// given u, we want compute v (L2–projection)
		// so we solve the following system:
		//
		//   cM . v = f, f := cfM . u
		//
		// for example, one can use CG (since mass matrix cM = cM^T and well–conditioned)
		// here we assemble these matrices given fine / coarse meshes and FE
		logger.beg("build pattern of mass matrices");
			  cM.generatePatternFrom(createDOFsConnectivityList(cMesh, FE));
			 cfM.generatePatternFrom(createDOFsConnectivityList(cMesh, fMesh, FE));
		logger.end();
		// master element
		Triangle2D master {
			{ {0., 0.}, {1., 0.}, {0., 1.} } 
		};
		// pieces of master element after refinement
		std::vector<Triangle2D> fElems {
			{ { { 0., 0. }, { .5, 0. }, { 0., .5 } } },
			{ { { 1., 0. }, { .5, .5 }, { .5, 0. } } },
			{ { { 0., 1. }, { 0., .5 }, { .5, .5 } } },
			{ { { .5, .5 }, { 0., .5 }, { .5, 0. } } }
		};
		// shapes
		auto masterShapes = FE.getShapesOf(master);
		// shapes for transformation T
		auto lagrangeMasterShapes = Triangle_P1_Lagrange::instance().getShapesOf(master);
		// current element and its middle nodes
		Triangle2D enodes;
		// jacobian of T
		double detGradT;
		// T, mapping from master element to physical ″
		VectorField2D T = [&](Node2D const & masterNode) {
			Node2D physicalNode { { 0., 0. } };
			for (LocalIndex i = 0; i < 3; ++i) {
				physicalNode[0] += enodes[i][0] * lagrangeMasterShapes[i](masterNode);
				physicalNode[1] += enodes[i][1] * lagrangeMasterShapes[i](masterNode);
			}
			return physicalNode;
		};
		// quagratures
		auto& qRule = MasterTriangleQuadratureRule::instance();
		LocalIndex deg = ceil(2. * FE.deg());
		logger.buf << "polynomial degree for Gaussian quadrature: " << deg;
		logger.log();
		// local matrices
		SymmetricMatrix<double> cLocalMassMatrixSkeleton(masterShapes.size());
		for (LocalIndex i = 0; i < masterShapes.size(); ++i)
			for (LocalIndex j = i; j < masterShapes.size(); ++j)
				cLocalMassMatrixSkeleton(i, j) = qRule.computeQuadrature([&](Node2D const & p) {
					return masterShapes[j](p) * masterShapes[i](p);
				}, deg);
		DenseMatrix<double> dummy(masterShapes.size());
		std::vector<decltype(dummy)> cfLocalMassMatrixSkeleton(fElems.size(), dummy);
		for (LocalIndex k = 0; k < fElems.size(); ++k) {
			enodes = fElems[k];
			for (LocalIndex i = 0; i < masterShapes.size(); ++i) 
				for (LocalIndex j = 0; j < masterShapes.size(); ++j)
					cfLocalMassMatrixSkeleton[k](i, j) += qRule.computeQuadrature([&](Node2D const & p) {
						return masterShapes[j](p) * masterShapes[i](T(p));
					}, deg);
		}
		logger.beg("assemble system and rhs mass matrices");
			for (Index ci = 0; ci < cMesh.numbOfElements(); ++ci) { // “ci” for coarse (element) index
				detGradT = 2. * area(cMesh.getElement(ci));
				// coarse dofs numn
				auto cDOFsNumn = FE.getDOFsNumeration(cMesh, ci);
				for (LocalIndex i = 0; i < cDOFsNumn.size(); ++i)
					for (LocalIndex j = i; j < cDOFsNumn.size(); ++j)
						cM(cDOFsNumn[i], cDOFsNumn[j]) += detGradT * cLocalMassMatrixSkeleton(i, j);
				// fine dofs numn for each fine element
				std::vector<decltype(cDOFsNumn)> fDOFsNumn;
				std::vector<Index> fElemsIndicies;
				fDOFsNumn.reserve(fElems.size());
				fElemsIndicies.reserve(fElems.size());
				for (Index fi : fMesh.getNeighborsIndicies(ci)) {
					fDOFsNumn.emplace_back(FE.getDOFsNumeration(fMesh, fi));
					fElemsIndicies.emplace_back(fi);
				}
				fDOFsNumn.emplace_back(FE.getDOFsNumeration(fMesh, ci));
				fElemsIndicies.emplace_back(ci);
				for (LocalIndex k = 0; k < fElems.size(); ++k) {
					detGradT = 2. * area(fMesh.getElement(fElemsIndicies[k]));
					for (LocalIndex i = 0; i < cDOFsNumn.size(); ++i)
						for (LocalIndex j = 0; j < cDOFsNumn.size(); ++j)
							cfM(cDOFsNumn[i], fDOFsNumn[k][j]) += detGradT * cfLocalMassMatrixSkeleton[k](i, j);
				}
			}
		logger.end();
		return {cM, cfM};
	}

}