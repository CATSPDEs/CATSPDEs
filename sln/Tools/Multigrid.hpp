#pragma once
#include "boost/tuple/tuple.hpp"
// logger
#include "SingletonLogger.hpp"
// mesh
#include "Triangulation.hpp"
// for prolongation matrices
#include "CSCMatrix.hpp"
// for the coarsest grid
#include "DenseMatrix.hpp"
// FEs
#include "AbstractFiniteElement.hpp"
// for grid transfer operators
#include "L2Projection.hpp"
#include "DivGradFEM.hpp" // mass matrix assembler
#include "ProjectionSolvers.hpp" // CG for easy mass matrix computations

/*
	Alexander Žilyakov, Jan 2017
*/

// helper types

template <typename TMatrix>
using Assembler = std::function<
	boost::tuple<TMatrix /* system matrix */, std::vector<double> /* rhs vector */>(Triangulation const &)
>;

template <typename TMatrix>
using Smoother = std::function<std::vector<double>(
	TMatrix& /* system matrix */, std::vector<double> const & /* rhs */, std::vector<double> const & /* initial guess */, Index /* mesh level */
)>;

enum class TransferType { canonical, L2 };

// functor for geometric multigrid routines

template <typename TMatrix>
class Multigrid {
	// data
	TransferType _transfer; // type of transfer operators (prolongations and restrictions)
	std::vector<double> _b; // fine rhs vector
	std::vector<TMatrix> _A; // system matrices
	DenseMatrix<double> _LU; // decomposition for coarse grid matrix
	std::vector<CSCMatrix<double>>  _P_Hh; // transfer operators matrices
	// for L2—projections transfer operators
	std::vector<CSlCMatrix<double>> _M_HH;
	std::vector<CSCMatrix<double>> _M_Hh;
	// helpers
	Index _numbOfDOFs(Index meshLevel) {
		return _A[meshLevel].getOrder();
	}
	std::vector<double> _prolongateTo(Index meshLevel, std::vector<double> const & u) {
		if (_transfer == TransferType::canonical)
			return _P_Hh[meshLevel - 1].t() * u;
		if (_transfer == TransferType::L2)
			return ProjectionSolvers::Krylov::CG(_M_HH[meshLevel], _M_Hh[meshLevel - 1].t() * u);
	}
	std::vector<double> _restrictFrom(Index meshLevel, std::vector<double> const & u) {
		if (_transfer == TransferType::canonical)
			return _P_Hh[meshLevel - 1] * u;
		if (_transfer == TransferType::L2)
			return ProjectionSolvers::Krylov::CG(_M_HH[meshLevel - 1], _M_Hh[meshLevel - 1] * u);
	}
public:
	Multigrid(
		TriangularScalarFiniteElement const & FE, // finite element
		Triangulation& Omega, // initial mesh
		Index numbOfMeshLevels, // numb of mesh levels
		Assembler<TMatrix> const & assembler, // system assembler
		TransferType transfer = TransferType::canonical  
	) : _transfer(transfer) {
		// logger
		auto& logger = SingletonLogger::instance();
		// reserve
		_A.reserve(numbOfMeshLevels + 1);
		if (transfer == TransferType::canonical) {
			_P_Hh.reserve(numbOfMeshLevels);
		}
		else if (transfer == TransferType::L2) {
			_M_HH.reserve(numbOfMeshLevels + 1);
			_M_Hh.reserve(numbOfMeshLevels);
		}
		logger.beg("mesh level 0");
			logger.beg("run assembler");
				auto system = assembler(Omega);
				_A.emplace_back(boost::get<0>(system));
			logger.end();
		logger.end();
		auto& fMesh = Omega; // synonim for fine mesh
		auto  cMesh = Omega; // copy for coarse mesh
		for (Index currentMeshLevel = 1; currentMeshLevel <= numbOfMeshLevels; ++currentMeshLevel) {
			logger.beg("mesh level " + std::to_string(currentMeshLevel));
				logger.beg("refine mesh");
					fMesh.refine();
				logger.end();
				logger.beg("run assembler");
					system = assembler(fMesh);
					_A.emplace_back(boost::get<0>(system));
				logger.end();
				if (transfer == TransferType::canonical) {
					logger.beg("assemble prolongation matrix");
						_P_Hh.emplace_back(FE.prolongate(cMesh, fMesh));
					logger.end();
				}
				else if (transfer == TransferType::L2) {
					logger.beg("assemble mass matrices for L2 projection transfer operators");
						auto mass = FEM::L2ProjectionAssembler(fMesh, cMesh, FE);
						_M_HH.emplace_back(static_cast<CSlCMatrix<double>>(boost::get<0>(mass)));
						_M_Hh.emplace_back(boost::get<1>(mass));
					logger.end();
				}
				cMesh = fMesh;
			logger.end();
		}
		if (transfer == TransferType::L2) {
			logger.beg("assemble the finest mass matrix for L2 projection transfer");
				// PDE, reaction term only
				ConvectionDiffusionEqn2D PDE {
					[](Node2D const &) { return 0.; },
					[](Node2D const &) { return 1.; },
					[](Node2D const &) { return 0.; }
				};
				ScalarBoundaryCondition2D RobinBC { {
						[](Node2D const &) { return 0.; },
						[](Node2D const &) { return 0.; } 
					},
					[](Node2D const & p) { return true; }
				}, DirichletBC {
					[](Node2D const &) { return 0.; }
				};
				_M_HH.emplace_back(boost::get<0>(
					FEM::DivGrad::assembleSystem(PDE, fMesh, RobinBC, DirichletBC, FE)
				));
			logger.end();
		}
		logger.beg("decompose coarse grid matrix");
			(_LU = _A.front()).decompose();
		logger.end();
		// fine rhs
		_b = boost::get<1>(system);
	}
	// MG iteration
	std::vector<double> operator()( // MG iteration for solving A.z = f
		Index meshLevel, // current mesh level
		std::vector<double> const & f, // rhs
		std::vector<double> const & z_0, // initial approximation to soln
		Smoother<TMatrix> const & smoother,
		Index gamma = 1 // numb of recursive calls (V-cycle by default)
	) {
		// logger
		auto& logger = SingletonLogger::instance();
		if (meshLevel) { // we are not at the coarsest grid
			auto& A = _A[meshLevel];
			logger.beg("mesh level " + std::to_string(meshLevel));
				logger.buf << "system size = " << A.getOrder();
				logger.log();
				logger.beg("smooth the residual (pre-smoothing)");
					//std::ofstream output("Mathematica/postprocessing/mg/vcycle.dat"/*, std::ios_base::app*/);
					//if (meshLevel == _A.size() - 1) {
					//	logger.wrn("export z_0");
					//	output << z_0 <<"\n";
					//}
					auto z = smoother(A, f, z_0, meshLevel);
					//if (meshLevel == _A.size() - 1) {
					//	logger.wrn("export smoothed z");
					//	output << z << "\n";
					//}
				logger.end();
				logger.beg("restrict the residual to mesh level " + std::to_string(meshLevel - 1));
					//std::ofstream output2("Mathematica/postprocessing/mg/residual.dat"/*, std::ios_base::app*/);
					//if (meshLevel == _A.size() - 1) {
					//	auto res = A * z - f;
					//	logger.buf << "  residual: " << res.size() << " " << norm(res) << "\n" << res << "\n"
					//	           << "r-residual: " << _restrictFrom(meshLevel, res).size() << " " << norm(_restrictFrom(meshLevel, res)) << "\n" << _restrictFrom(meshLevel, res);
					//	logger.log();
					//	output2 << res << '\n' << _restrictFrom(meshLevel, res);
					//}
					auto d = _restrictFrom(meshLevel, A * z - f);
				logger.end();
				logger.log("go to coarser grid");
					std::vector<double> e(_numbOfDOFs(meshLevel - 1)); // initial guess for the error on the coarser mesh level
					for (Index i = 0; i < gamma; ++i)
						e = operator()(meshLevel - 1, d, e, smoother, gamma);
				logger.beg("correct from the coarse grid");
					//std::ofstream output1("Mathematica/postprocessing/mg/correction.dat"/*, std::ios_base::app*/);
					//if (meshLevel == _A.size() - 1) {
					//	logger.buf << "  correction: " << e.size() << " " << norm(e) << "\n" << e << "\n"
					//	           << "p-correction: " << _prolongateTo(meshLevel, e).size() << " " << norm(_prolongateTo(meshLevel, e)) << "\n" << _prolongateTo(meshLevel, e);
					//	logger.log();
					//	output1 << e << '\n' << _prolongateTo(meshLevel, e);
					//}
					z -= _prolongateTo(meshLevel, e);
					//if (meshLevel == _A.size() - 1) {
					//	logger.wrn("export corrected z");
					//	output << z <<"\n";
					//}
				logger.end();
				logger.beg("smooth the residual (post-smoothing)");
					z = smoother(A, f, z, meshLevel);
					//if (meshLevel == _A.size() - 1) {
					//	logger.wrn("export post-smoothed z");
					//	output << z <<"\n";
					//}
				logger.end();
			logger.end();
			return z;
		}
		// we are at the coarsest grid
		logger.beg("mesh level 0");
			DenseMatrix<double> A { _A.front() };
			logger.buf << "system size = " << A.getOrder();
			logger.log();
			logger.beg("compute exact soln");
				// auto z = A.GaussElimination(f);
				auto z = _LU.backSubst(_LU.forwSubst(f, 0.));
				logger.buf << "residual norm = " << norm(_A.front() * z - f);
				logger.log();
			logger.end();
		logger.end();
		return z;
	}
	// get system matrices
	auto& A() { return _A.back(); }
	auto& A(Index meshLevel) { return _A[meshLevel]; }
	// get fine rhs vector
	auto& b() { return _b; }
	// for analysis
	auto& P_Hh(Index fromMeshLevel) { return _P_Hh[fromMeshLevel]; }
};