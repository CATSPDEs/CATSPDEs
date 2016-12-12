#include "FEM.hpp"
#include "ProjectionSolvers.hpp" // conjugate gradients
#include "SingletonLogger.hpp"

using std::string;
using std::to_string;
using namespace FEM::DivGrad;

// logger
SingletonLogger& logger = SingletonLogger::instance();

int main() {
	auto eps = 10e-8;
	string oPath("Mathematica/dummy/res/");
	try {
		logger.beg("set up PDE, BCs, and import initial mesh");
			DiffusionReactionEqn2D PDE(
				[](Node2D const & p) { return 1.; },
				[](Node2D const & p) { return 0.; },
				[](Node2D const & p) { return 2.*(-1. + p[0])*p[0] + 2.*(-1. + p[1])*p[1]; }
			);
			Triangulation Omega;
			Omega.import("mesh.ntn");
			DirichletScalarCondition2D DirichletBC(
				[](Node2D const & p) { return 1.; }
			);
		logger.end();
		Index numbOfMeshLevels;
		logger.inp("numb of mesh levels", numbOfMeshLevels);
		auto method = logger.opt("choose solving technique", {
			"MG",
			"PCG w/ MG as an inner iteration",
			"CG",
			"Jacobi",
			"Gauss-Seidel",
			"Collect data for different mesh levels (system size, nonzeros, etc.)",
			"MG iteration complexity"
		});
		if (method == 0) {
			using namespace Multigrid;
			logger.beg("setup multigrid data");
				auto system = linearLagrangeSetter(PDE, DirichletBC, Omega, numbOfMeshLevels);
				auto& A = boost::get<0>(system);
				auto& b = boost::get<1>(system);
				Index maxNumbOfIterations;
				logger.inp("max numb of iterations", maxNumbOfIterations);
				Index smootherType = logger.opt("smoothing technique", { "Jacobi", "Gauss-Seidel" });
				double omega;
				logger.inp("relaxation parameter", omega);
				Index nu;
				logger.inp("numb of smoothing iterations", nu);
				if (smootherType == 0) // Jacobi
					Multigrid::smoother = [&](
							AbstractMultipliableMatrix<double>& A, // system matrix
							std::vector<double> const & b, // rhs
							std::vector<double> const & x_0 // initial guess
						) {
							return boost::get<0>(ProjectionSolvers::Smoothers::Jacobi(A, b, x_0, omega, nu));
						};
				else // Gauss-Seidel
					Multigrid::smoother = [&](
							SymmetricCSlCMatrix<double>& A, // system matrix
							std::vector<double> const & b, // rhs
							std::vector<double> const & x_0 // initial guess
						) {
							return boost::get<0>(ProjectionSolvers::Smoothers::GaussSeidel(A, b, x_0, omega, nu));
						};
				gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
				++gamma;
			logger.end();
			logger.beg("solve w/ MG");
				std::vector<double> xi(A.getOrder(), 0.), r;
				for (Index i = 0; i < maxNumbOfIterations; ++i) {
					r.emplace_back(norm(b - A * xi));
					logger.buf << "|| r_" << i << " || = " << std::scientific << r.back();
					logger.log();
					if (r.back() < eps) break;
					logger.mute = true;
					xi = iteration(numbOfMeshLevels, xi, b);
					logger.mute = false;
				}
			logger.end();
			logger.beg("export residuals norms history");
				//Omega.export(oPath + "mesh.nt");
				export(xi, oPath + "xi.dat");
				export(r, oPath + "MG.dat");
			logger.end();
		} 
		else if (method == 1) {
			using namespace Multigrid;
			logger.beg("setup multigrid data");
				auto system = linearLagrangeSetter(PDE, DirichletBC, Omega, numbOfMeshLevels);
				auto& A = boost::get<0>(system);
				auto& b = boost::get<1>(system);
				Index numbOfInnerIterations;
				logger.inp("numb of iterations for inner solver", numbOfInnerIterations);
				Index smootherType = logger.opt("smoothing technique", { "Jacobi", "Gauss-Seidel" });
				double omega;
				logger.inp("relaxation parameter", omega);
				Index nu;
				logger.inp("numb of smoothing iterations", nu);
				if (smootherType == 0) // Jacobi
					Multigrid::smoother = [&](
							AbstractMultipliableMatrix<double>& A, // system matrix
							std::vector<double> const & b, // rhs
							std::vector<double> const & x_0 // initial guess
						) {
							return boost::get<0>(ProjectionSolvers::Smoothers::Jacobi(A, b, x_0, omega, nu));
						};
				else // Gauss-Seidel
					Multigrid::smoother = [&](
							SymmetricCSlCMatrix<double>& A, // system matrix
							std::vector<double> const & b, // rhs
							std::vector<double> const & x_0 // initial guess
						) {
							return boost::get<0>(ProjectionSolvers::Smoothers::GaussSeidel(A, b, x_0, omega, nu));
						};
				gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
				++gamma;
			logger.end();
			logger.beg("solve w/ PCG");
				auto result = ProjectionSolvers::Krylov::PCG(A, b, [&](std::vector<double> const & x) {
					std::vector<double> y(A.getOrder(), 0.);
					logger.mute = true;
					for (Index i = 0; i < numbOfInnerIterations; ++i)
						y = iteration(numbOfMeshLevels, y, x);
					logger.mute = false;
					return y;
				}, boost::none, eps);
				auto xi = boost::get<0>(result);
				auto r = boost::get<1>(result);
			logger.end();
			logger.beg("export residuals norms history");
				export(xi, oPath + "xi.dat");
				export(r, oPath + "PCG.dat");
			logger.end();
		}
		else if (method == 2) {
			logger.beg("refine mesh");
				Omega.refine(numbOfMeshLevels);
			logger.end();
			logger.beg("assemble system");
				auto system = linearLagrangeAssembler(PDE, Omega, DirichletBC);
				auto& A = boost::get<0>(system);
				auto& b = boost::get<1>(system);
			logger.end();
			logger.beg("solve w/ CG");
				auto result = ProjectionSolvers::Krylov::CG(A, b, boost::none, eps);
				auto xi = boost::get<0>(result);
				auto r = boost::get<1>(result);
				auto z = norm(b - A * xi);
				if (z >= eps) logger.wrn("CG failed");
				logger.buf << "|| r_" << r.size() - 1 << " || = " << std::scientific << z;
				logger.log();
			logger.end();
			logger.beg("export residuals norms history");
				export(xi, oPath + "xi.dat");
				export(r, oPath + "CG.dat");
			logger.end();
		}
		else if (method == 3) {
			logger.beg("refine mesh");
				Omega.refine(numbOfMeshLevels);
			logger.end();
			logger.beg("assemble system");
				auto system = linearLagrangeAssembler(PDE, Omega, DirichletBC);
				auto& A = boost::get<0>(system);
				auto& b = boost::get<1>(system);
			logger.end();
			logger.beg("setup solver data");
				Index maxNumbOfIterations;
				logger.inp("max numb of iterations", maxNumbOfIterations);
				double omega;
				logger.inp("relaxation parameter", omega);
			logger.end();
			logger.beg("solve w/ Jacobi");
				auto result = ProjectionSolvers::Smoothers::Jacobi(A, b, boost::none, omega, maxNumbOfIterations);
				auto r = boost::get<1>(result);
			logger.end();
			logger.beg("export residuals norms history");
				export(r, oPath + "Jacobi.dat");
			logger.end();
		}
		else if (method == 4) {
			logger.beg("refine mesh");
				Omega.refine(numbOfMeshLevels);
			logger.end();
			logger.beg("assemble system");
				auto system = linearLagrangeAssembler(PDE, Omega, DirichletBC);
				auto& A = boost::get<0>(system);
				auto& b = boost::get<1>(system);
			logger.end();
			logger.beg("setup solver data");
				Index maxNumbOfIterations;
				logger.inp("max numb of iterations", maxNumbOfIterations);
				double omega;
				logger.inp("relaxation parameter in (0., 2.)", omega);
			logger.end();
			logger.beg("solve w/ Gauss-Seidel");
				auto result = ProjectionSolvers::Smoothers::GaussSeidel(A, b, boost::none, omega, maxNumbOfIterations);
				auto r = boost::get<1>(result);
			logger.end();
			logger.beg("export residuals norms history");
				export(r, oPath + "GS.dat");
			logger.end();
		}
		else if (method == 5) {
			std::vector<Index> meshLevels, systemSizes, nnzeros;
			for (Index i = 1; i <= numbOfMeshLevels; ++i) {
				logger.beg("mesh level " + to_string(i));
					logger.beg("refine mesh");
						Omega.refine();
					logger.end();
					logger.beg("assemble system");
						auto system = linearLagrangeAssembler(PDE, Omega, DirichletBC);
						auto& A = boost::get<0>(system);
					logger.end();
					meshLevels.push_back(i);
					systemSizes.push_back(A.getOrder());
					nnzeros.push_back(A.nnz());
				logger.end();
			}
			logger.beg("export data");
				export(meshLevels, oPath + "meshLevels.dat");
				export(systemSizes, oPath + "systemSizes.dat");
				export(nnzeros, oPath + "nnzeros.dat");
			logger.end();
		}
		else if (method == 6) {
			using namespace Multigrid;
			logger.beg("setup multigrid data");
				auto system = linearLagrangeSetter(PDE, DirichletBC, Omega, numbOfMeshLevels);
				auto& A = boost::get<0>(system);
				auto& b = boost::get<1>(system);
				Index smootherType = logger.opt("smoothing technique", { "Jacobi", "Gauss-Seidel" });
				Index nu;
				logger.inp("numb of smoothing iterations", nu);
				if (smootherType == 0) // Jacobi
					Multigrid::smoother = [&](
							AbstractMultipliableMatrix<double>& A, // system matrix
							std::vector<double> const & b, // rhs
							std::vector<double> const & x_0 // initial guess
						) {
							return boost::get<0>(ProjectionSolvers::Smoothers::Jacobi(A, b, x_0, .1, nu));
						};
				else // Gauss-Seidel
					Multigrid::smoother = [&](
							SymmetricCSlCMatrix<double>& A, // system matrix
							std::vector<double> const & b, // rhs
							std::vector<double> const & x_0 // initial guess
						) {
							return boost::get<0>(ProjectionSolvers::Smoothers::GaussSeidel(A, b, x_0, 1., nu));
						};
				gamma = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" });
				++gamma;
			logger.end();
			logger.beg("try 10 MG iterations");
				logger.mute = true;
				std::vector<double> xi(A.getOrder(), 0.), r;
				for (Index i = 0; i < 10; ++i)
					xi = iteration(numbOfMeshLevels, xi, b);
				logger.mute = false;
			logger.end();
		}
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}