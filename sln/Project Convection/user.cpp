#include <string>
#include <numeric>
// finite elements to use
#include "Triangle_P1_Lagrange.hpp"
// interpolant for diffusion
#include "FEInterpolant.hpp"
// assembler
#include "DivGradFEM.hpp"
// solvers
#include "ProjectionSolvers.hpp"
// pi
#include "constants.hpp"

using std::string;
using std::vector;
using boost::get;
using namespace FEM::DivGrad;
using namespace ProjectionSolvers;

int main() {
	string iPath("Mathematica/preprocessing/"), oPath("Mathematica/postprocessing/");
	auto& logger = SingletonLogger::instance();
	try {
		logger.beg("import mesh");
			Triangulation Omega;
			Omega.import(iPath + "mesh.nt");
			Omega.computeNeighbors().enumerateRibs();
		logger.end();
		logger.beg("set up eqn");
			Quadrilateral2D const water {
				{ {0., 0.}, {.2, 0.}, {.2, .2}, {0., .2} }
			};
			auto const c = [&](Node2D const & p) {
				return nodeInElement(water, p) ? 4180. : 540.;
			};
			auto const rho = [&](Node2D const & p) {
				return nodeInElement(water, p) ? 10. : 7000.;
			};
			auto const lambda = [&](Node2D const & p) {
				return nodeInElement(water, p) ? .6 : 47.2;
			};
			double u_air, u_ini;
			logger.inp("set initial temperature of water", u_ini);
			logger.inp("set outside temperature", u_air);
			// bc
			double heatFlux;
			logger.inp("set heat flux [Watt]", heatFlux);
			heatFlux /= (PI * .21 * .21); // per area
			double const R = 5.7; // for Robin coef (iron-air)
			auto approxEqual = [](double x, double y) { return fabs(x - y) < 1e-12; };
			ScalarField2D NeumannValue = [&](Node2D const & p) {
				if (approxEqual(p[0], 0.)) return 0.;
				if (approxEqual(p[1], -.01)) return heatFlux;
				return R * u_air;
			};
			ScalarField2D RobinCoefficient = [&](Node2D const & p) {
				if (approxEqual(p[0], 0.) || approxEqual(p[1], -.01)) return 0.;
				return R;
			};
			ScalarBoundaryCondition2D
				RobinBC {
					{ RobinCoefficient, NeumannValue }, // n . (a ∇u) + R u = N
					[](Node2D const &) { return true; }
				},
				DirichletBC { [](Node2D const &) { return 0.;} };
			// wind
			double const wind_const = 4.37e+8;
			double wind_mean;
			logger.inp("set mean wind velocity", wind_mean);
			wind_mean *= wind_const;
			auto const wind = [&](Node2D const & p) -> Node2D {
				if (nodeInElement(water, p)) {
					Node2D res {
						-2. * (-.2 + p[0]) * p[0] * p[0] * (.02 + p[0] * (-.1 + p[1]) - .2 * p[1]) * (-.2 + p[1]) * p[1],
						2. * (-.2 + p[0]) * p[0] * (.02 + p[0] * (-.2 + p[1]) - .1 * p[1]) * (-.2 + p[1]) * p[1] * p[1]
					};
					return wind_mean * res;
				}
				return {0., 0.};
			};
			// initial data
			auto& FE = Triangle_P1_Lagrange::instance();
			vector<vector<double>> solns;
			solns.push_back(vector<double>(FE.numbOfDOFs(Omega), u_ini));
		logger.end();
		logger.beg("start time steps");
			double delta_t;
			logger.inp("set time step", delta_t);
			Index n;
			logger.inp("set max numb of time steps", n);
			// iterations logger
			Index i_log;
			logger.inp("log every nth iteration of the solver, n", i_log);
			// max numb of iters
			Index maxNumbOfIterations;
			logger.inp("max numb of iterations", maxNumbOfIterations);
			// tolerance
			double eps;
			logger.inp("set eps for solver", eps);
			// save dofs indicies in water
			std::set<Index> waterIndicies;
			auto centerWater = [](Node2D const & p) {
				return (p[0] - .1) * (p[0] - .1) + (p[1] - .1) * (p[1] - .1) <= .05 * .05;
			};
			for (Index t = 0; t < Omega.numbOfElements(); ++t) {
				auto nodes = FE.getDOFsNodes(Omega, t);
				auto indicies = FE.getDOFsNumeration(Omega, t);
				for (LocalIndex k = 0; k < nodes.size(); ++k)
					if (centerWater(nodes[k])) waterIndicies.insert(indicies[k]);
			}
			for (Index i = 1; i <= n; ++i) {
				logger.buf
					<< "time frame #:   " << i << '\n'
					<< "time [seconds]: " << i * delta_t;
				logger.log();
				TriangularScalarFEInterpolant u_prev { solns[i - 1], FE, Omega };
				Index activeElementIndex;
				ScalarField2D 
					diffusion = [&](Node2D const & p) { 
						return delta_t * lambda(p);
					},
					reaction = [&](Node2D const & p) { 
						return c(p) * rho(p); 
					},
					force = [&](Node2D const & p) {
						return c(p) * rho(p) * u_prev(p, activeElementIndex);
					};
					VectorField2D convection = [&](Node2D const & p) -> Node2D {
						//return { 0., 0. };
						return delta_t * c(p) * rho(p) * wind(p);
					};
				DiffusionReactionEqn2D PDE { diffusion, reaction, force, convection };
				logger.beg("assemble system");
					auto system = assembleSystem(PDE, Omega, RobinBC, DirichletBC, FE, activeElementIndex);
					auto& A = get<0>(system);
					auto& b = get<1>(system);
				logger.end();
				logger.beg("build ILU(0)");
					auto LDLT = A;
					LDLT.decompose();
					auto P = [&](vector<double> const & x) {
						return LDLT.backSubst(LDLT.diagSubst(LDLT.forwSubst(x, 0.)), 0.);
					};
				logger.end();
				logger.beg("solve linear problem");
					solns.push_back(
						Krylov::PBiCGStab(P, A, b, solns[i - 1], maxNumbOfIterations, eps, StoppingCriterion::relative, i_log)
					);
				logger.end();
				// mean temperature
				double meanTemp = 0.;
				for (Index k : waterIndicies) meanTemp += solns.back()[k];
				meanTemp /= waterIndicies.size();
				logger.buf << "mean temperature: " << meanTemp;
				logger.log();
				if (meanTemp >= 100.) {
					logger.log("100 degrees!");
					break;
				}
				// if (accumulate(solns.back().begin(), solns.back().end(), 0.) / solns.back().size() >= 100.) break;
			}
			logger.beg("export solns");
				export(solns, oPath + "solns_wind.dat");
			logger.end();
		logger.end();
		logger.exp("stdin.txt");
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}