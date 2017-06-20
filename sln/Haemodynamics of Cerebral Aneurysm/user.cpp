#include <memory>
// Taylor–Hood
#include "Triangle_P2_Lagrange.hpp"
#include "Triangle_P1_Lagrange.hpp"
// MINI
#include "Triangle_Pt3_LagrangeBubble.hpp"
#include "Triangle_P1_Lagrange.hpp"
// assemblers
#include "MixedFEM.hpp"
#include "DivGradFEM.hpp"
// saddle point system
#include "SymmetricBlockMatrix.hpp" 
// solvers
#include "ProjectionSolvers.hpp"
// precond
#include "Multigrid.hpp"
// interpolant for wind field
#include "FEInterpolant.hpp"
// interpolant for inflow amplitude
#include "SegmentMesh.hpp"
#include "Segment_P3_Hermite.hpp"

using std::vector;
using std::string;
using std::shared_ptr;
using boost::get;
using namespace ProjectionSolvers;
using namespace FEM;

int main() {
	string iPath("Mathematica/preprocessing/"), oPath("Mathematica/postprocessing/");
	auto& logger = SingletonLogger::instance();
	try {
		// current time and time step
		double t = 0., delta_t;
		Index final_frame_index;
		logger.inp("set time step", delta_t);
		logger.inp("final time frame number", final_frame_index);
		logger.beg("set up mesh");
			vector<string> meshType { "mesh_fine.ntn", "mesh_coarse.ntn" };
			auto meshTypeIndex = logger.opt("mesh type", meshType);
			Triangulation Omega_0;
			Omega_0.import(iPath + meshType[meshTypeIndex]);
			Omega_0.enumerateRibs();
			auto Omega = Omega_0;
			struct MeshData {
				Index numb_of_levels;
			} meshData;
			logger.inp("numb of mesh levels", meshData.numb_of_levels);
		logger.end();
		logger.beg("set up units");
			// blood parameters
			struct BloodParameters {
				double mu, rho; // pascals-seconds, kgs per meter cubed
				double v_min, v_max; // meters
			} const bloodParameters {
				3.3e-03, // dynamic viscocity
				1060, // density
				.7, 1.5 // systolic and diastolic velocities
			};
			// geometry
			struct VesselsGeometry {
				double length, radius;
			} const vesselsGeometry { 16., .002 };
			// units
			struct Units {
				double length, velocity, time; // meters, meters per second, seconds
			} const units { 
				vesselsGeometry.radius, // radius of the aurtery
				bloodParameters.v_max, // systolic velocity
				.875 // cardiac cycle
			};
			// Reynolds number

			double const Re = units.velocity * units.length * bloodParameters.rho / bloodParameters.mu;
			//double const Re = 100;
			
			// Strouhal number
			double const St = units.length / (units.time * units.velocity);
			logger.buf
				<< "Reynolds number: " << Re << '\n'
				<< "Strouhal number: " << St;
			logger.log();
		logger.end();
		logger.beg("set up PDE and BCs");
			auto approxEqual = [](double x, double y) { return fabs(x - y) < 1e-12; };
			// Stokes problem for (u_0, p_0)
			OseenProblem2D Oseen {
				0.,
				[](Node2D const &) -> Node2D { return { 0., 0. }; },
				1. / Re,
				[](Node2D const &) -> Node2D { return { 0., 0. }; },
				[](Node2D const &) { return 0.; }
			};
			// natural BCs for outflow boundary
			auto GammaOut = [&](Node2D const & p) { return approxEqual(p[0], vesselsGeometry.length / 2.); };
			auto h = [&](Node2D const &) -> Node2D { return { 0., 0.}; };
			VectorBoundaryCondition2D NeumannBC { h, GammaOut };
			// essential BCs for inflow boundary and no-slip boundary
			auto GammaIn = [&](Node2D const & p) { return approxEqual(p[0], -vesselsGeometry.length / 2.); };
			// systolic / diastolic amplitude
			SegmentMesh hermiteNodes {
				{ 0., units.time / 7., units.time }
			};
			SegmentFEInterpolant interp { {
					bloodParameters.v_min,
					bloodParameters.v_max,
					bloodParameters.v_min,
					0., 0., 0.
				}, Segment_P3_Hermite::instance(), hermiteNodes 
			};
			SmartScalarField1D ampl = [&](double const & t) {
				// s := units.time * t, [s] = seconds
				return interp(fmod(units.time * t, units.time)) / bloodParameters.v_max;
			};
			auto g = [&](Node2D const & p) -> Node2D { 
				if (GammaIn(p)) return { ampl(t) * (1 - p[1] * p[1]), 0. };
				return { 0., 0. };
			};
			VectorBoundaryCondition2D DirichletBC { g };
		logger.end();
		auto FEPairIndex = logger.opt("choose FE pair", { "Taylor-Hood", "MINI" });
		// Taylor—Hood family
		TriangularScalarFiniteElement
			*velocityFE = &Triangle_P2_Lagrange::instance(),
			*pressureFE = &Triangle_P1_Lagrange::instance();
		// MINI element
		if (FEPairIndex == 1) {
			velocityFE = &Triangle_Pt3_LagrangeBubble::instance();
			pressureFE = &Triangle_P1_Lagrange::instance();
		};
		logger.beg("set solver options");
			struct SolverOptions {
				Index outer_iters_max, inner_iters, smoothing_iters, precond_index, i_log, cycle, smoother_index;
				double eps;
				StoppingCriterion stop;
			} opt;
			opt.precond_index = logger.opt("choose precond for BiCGStab", { "I", "P_BD", "P_BT" });
			logger.inp("max numb of iterations for outer solver", opt.outer_iters_max);
			logger.inp("log every nth iteration, n", opt.i_log);
			logger.inp("set eps for solver", opt.eps);
			opt.stop = (StoppingCriterion)logger.opt("stopping criterion", { "absolute", "relative" });
			logger.inp("numb of iterations for inner solver", opt.inner_iters);
			opt.smoother_index = logger.opt("choose smoother", { "forwSOR", "SSOR", "ILU(0)" });
			logger.inp("numb of smoothing iters", opt.smoothing_iters);
			opt.cycle = logger.opt("set recursive calls type", { "V-cycle", "W-cycle" }) + 1;
		logger.end();
		logger.beg("build preconditioner and assemble system");
			Preconditioner P;
			struct {
				CSlCMatrix<double> A11;
				CSCMatrix<double> B1, B2;
				vector<double> b;
			} system;
			Index n, m;
			shared_ptr<Multigrid<CSlCMatrix<double>>> MG;
			SymmetricCSlCMatrix<double> pressureMassMatrix;
			Smoother<CSlCMatrix<double>> smoother;
			vector<CSlCMatrix<double>> ILU0;
			ILU0.reserve(meshData.numb_of_levels);
			if (opt.precond_index == 0) {
				P = [](vector<double> const & x) { return x; }; // do nothing
				logger.beg("refine mesh");
					Omega.refine(meshData.numb_of_levels);
				logger.end();
				auto res = Mixed::assembleSystem(Oseen, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE);
				system.A11 = get<0>(res);
				system.B1  = get<1>(res);
				system.B2  = get<2>(res);
				system.b   = get<3>(res);
			}
			else {
				logger.beg("(1) MG for Laplace block");
					MG = std::make_shared<Multigrid<CSlCMatrix<double>>>(
						Multigrid<CSlCMatrix<double>> {
							*velocityFE, Omega, meshData.numb_of_levels,
							[&](Triangulation const & Omega) {
								auto res = Mixed::assembleSystem(Oseen, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE);
								system.A11 = get<0>(res);
								system.B1  = get<1>(res);
								system.B2  = get<2>(res);
								system.b   = get<3>(res);
								return system.A11; // return Laplace block
							}
						}
					);
				logger.end();
				logger.beg("(2) pressure mass matrix for Schur complement");
					DiffusionReactionEqn2D ReactionEqn {
						[](Node2D const &) { return 0.; },
						[](Node2D const &) { return 1.; },
						[](Node2D const &) { return 0.; }
					};
					ScalarBoundaryCondition2D rBC { {
								[](Node2D const &) { return 0.; },
								[](Node2D const &) { return 0.; }
							},
							[](Node2D const & p) { return true; }
						}, dBC { [](Node2D const &) { return 0.; } };
					pressureMassMatrix = get<0>(DivGrad::assembleSystem(ReactionEqn, Omega, rBC, dBC, *pressureFE));
				logger.end();
				logger.beg("define smoother");
					if (opt.smoother_index == 0)
						smoother = [&](CSlCMatrix<double>& A, vector<double> const & b, vector<double> const & x_0, Index meshLevel) {
							return Smoothers::forwSOR(A, b, x_0, 1., opt.smoothing_iters, 0., StoppingCriterion::absolute, 0);
						};
					else if (opt.smoother_index == 1)
						smoother = [&](CSlCMatrix<double>& A, vector<double> const & b, vector<double> const & x_0, Index meshLevel) {
							return Smoothers::SSOR(A, b, x_0, 1., opt.smoothing_iters, 0., StoppingCriterion::absolute, 0);
						};
					else if (opt.smoother_index == 2) {
						logger.beg("compute ILU0's");
							for (Index i = 0; i < meshData.numb_of_levels; ++i) {
								auto ILU0i = MG->A(i + 1);
								ILU0.emplace_back(ILU0i.decompose());
							}
						logger.end();
						smoother = [&](CSlCMatrix<double>& A, vector<double> const & b, vector<double> const & x_0, Index meshLevel) {
							auto& d = ILU0[meshLevel - 1];
							return d.backSubst(d.diagSubst(d.forwSubst(b, 0.)), 0.);
						};
					}
				logger.end();
				if (opt.precond_index == 1) // BD
					P = [&](vector<double> const & x) {
						vector<double>
							x1(x.begin(), x.begin() + n),
							x2(x.begin() + n, x.begin() + 2 * n),
							x3(x.begin() + 2 * n, x.end()),
							y1(n), y2(n);
						// (1) Laplace block
						logger.mute = true;
						for (Index i = 0; i < opt.inner_iters; ++i) {
							y1 = (*MG)(meshData.numb_of_levels, x1, y1, smoother, opt.cycle);
							y2 = (*MG)(meshData.numb_of_levels, x2, y2, smoother, opt.cycle);
						}
						logger.mute = false;
						// (2) Schur complement
						auto y3 = pressureMassMatrix.diagSubst(x3);
						y3 *= -1.;
						// final vector
						vector<double> y;
						y.reserve(2 * n + m);
						y.insert(y.end(), y1.begin(), y1.end());
						y.insert(y.end(), y2.begin(), y2.end());
						y.insert(y.end(), y3.begin(), y3.end());
						return y;
					};
				else if (opt.precond_index == 2) // BT
					P = [&](vector<double> const & x) {
						vector<double>
							x1(x.begin(), x.begin() + n),
							x2(x.begin() + n, x.begin() + 2 * n),
							x3(x.begin() + 2 * n, x.end()),
							y1(n), y2(n);
						// (1) Schur complement
						auto y3 = pressureMassMatrix.diagSubst(x3);
						y3 *= -1.;
						// (1) Laplace block
						auto f1 = x1 - system.B1.t() * y3;
						auto f2 = x2 - system.B2.t() * y3;
						logger.mute = true;
						for (Index i = 0; i < opt.inner_iters; ++i) {
							y1 = (*MG)(meshData.numb_of_levels, f1, y1, smoother, opt.cycle);
							y2 = (*MG)(meshData.numb_of_levels, f2, y2, smoother, opt.cycle);
						}
						logger.mute = false;
						// final vector
						vector<double> y;
						y.reserve(2 * n + m);
						y.insert(y.end(), y1.begin(), y1.end());
						y.insert(y.end(), y2.begin(), y2.end());
						y.insert(y.end(), y3.begin(), y3.end());
						return y;
					};
			}
			n = system.A11.getOrder();
			m = system.B1.numbOfRows();
			SymmetricBlockMatrix<double> A {
				{ &system.A11, &system.A11, nullptr }, // diag
				{ // lval
					nullptr,
					&system.B1, &system.B2
				}, BlockSymmetryType::antisymmetric
			};
			for (Index i = 2 * n; i < system.b.size(); ++i) system.b[i] = -system.b[i];
		logger.end();
		logger.beg("Stokes problem for (u_0, p_0)");
			logger.beg("solve");
				auto x_0 = Krylov::BiCGStab(A, system.b, boost::none, 10, 1., StoppingCriterion::absolute, 1);
				x_0 = Krylov::PBiCGStab(
					P, A, system.b, x_0,
					opt.outer_iters_max, opt.eps, opt.stop, opt.i_log
				);
			logger.end();
			logger.beg("export soln vector");
				export(x_0, oPath + "x_0.dat");
			logger.end();
			auto exportMatrix = logger.yes("export matrix blocks");
			if (exportMatrix) {
				logger.beg("export blocks of matrix and rhs vector");
					static_cast<CSCMatrix<double>>(system.A11).exportHarwellBoeing(oPath + "system/A11_0.rua");
					system.B1.exportHarwellBoeing(oPath + "system/B1_0.rra");
					system.B2.exportHarwellBoeing(oPath + "system/B2_0.rra");
					export(system.b, oPath + "system/b_0.dat");
				logger.end();
			}
		logger.end();
		if (final_frame_index) {
			logger.beg("Oseen Problem for (u_1, p_1)");
				logger.beg("update parameters");
					t = delta_t;
					Oseen.massCoef() = St / delta_t;
					logger.buf << "alpha = " << Oseen.massCoef();
					logger.log();
					// wind
					Index activeElementIndex;
					TriangularScalarFEInterpolant 
						u_0_x_interp { std::vector<double>(x_0.begin(), x_0.begin() + n), *velocityFE, Omega },
						u_0_y_interp { std::vector<double>(x_0.begin() + n, x_0.begin() + 2 * n), *velocityFE, Omega };
					Oseen.windField() = [&](Node2D const & p) -> Node2D {
						return { u_0_x_interp(p, activeElementIndex), u_0_y_interp(p, activeElementIndex) };
					};
					// force
					Oseen.forceTerm() = [&](Node2D const & p) -> Node2D {
						return { St * u_0_x_interp(p, activeElementIndex) / delta_t, St * u_0_y_interp(p, activeElementIndex) / delta_t };
					};
				logger.end();
				if (opt.precond_index == 0) {
					auto res = Mixed::assembleSystem(Oseen, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE, activeElementIndex);
					system.A11 = get<0>(res);
					system.B1  = get<1>(res);
					system.B2  = get<2>(res);
					system.b   = get<3>(res);
				}
				else {
					logger.beg("REDEFINE MG for Laplace block");
						Omega = Omega_0;
						MG = std::make_shared<Multigrid<CSlCMatrix<double>>>(
							Multigrid<CSlCMatrix<double>> {
								*velocityFE, Omega, meshData.numb_of_levels,
								[&](Triangulation const & Omega) {
									auto res = Mixed::assembleSystem(Oseen, Omega, NeumannBC, DirichletBC, *velocityFE, *pressureFE, activeElementIndex);
									system.A11 = get<0>(res);
									system.B1  = get<1>(res);
									system.B2  = get<2>(res);
									system.b   = get<3>(res);
									return system.A11; // return Laplace block
								}
							}
						);
					logger.end();
					if (opt.smoother_index == 2) {
						logger.beg("REDEFINE smoother");
							ILU0.resize(0);
							ILU0.reserve(meshData.numb_of_levels);
							logger.beg("compute ILU0's");
								for (Index i = 0; i < meshData.numb_of_levels; ++i) {
									auto ILU0i = MG->A(i + 1);
									ILU0.emplace_back(ILU0i.decompose());
								}
							logger.end();
							smoother = [&](CSlCMatrix<double>& A, vector<double> const & b, vector<double> const & x_0, Index meshLevel) {
								auto& d = ILU0[meshLevel - 1];
								return d.backSubst(d.diagSubst(d.forwSubst(b, 0.)), 0.);
							};
						logger.end();
					}
				}
				for (Index i = 2 * n; i < system.b.size(); ++i) system.b[i] = -system.b[i];
				logger.beg("solve");
					auto x_1 = Krylov::BiCGStab(A, system.b, x_0, 10, 1., StoppingCriterion::absolute, 1);
					x_1 = Krylov::PBiCGStab(
						P, A, system.b, x_1, 
						opt.outer_iters_max, opt.eps, opt.stop, opt.i_log
					);
				logger.end();
				logger.beg("export soln vector");
					export(x_1, oPath + "x_1.dat");
				logger.end();
				if (exportMatrix) {
					logger.beg("export blocks of matrix and rhs vector");
						static_cast<CSCMatrix<double>>(system.A11).exportHarwellBoeing(oPath + "system/A11_1.rua");
						system.B1.exportHarwellBoeing(oPath + "system/B1_1.rra");
						system.B2.exportHarwellBoeing(oPath + "system/B2_1.rra");
						export(system.b, oPath + "system/b_1.dat");
					logger.end();
				}
			logger.end();
		}
		logger.beg("export mesh");
			Omega.export(oPath + "mesh.ntnr", { { "format", "NTNR" } });
		logger.end();
		logger.exp("stdin.txt");
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}