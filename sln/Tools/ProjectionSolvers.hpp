#pragma once
#include <boost/tuple/tuple.hpp>
#include <boost/optional/optional.hpp>
#include "AbstractPreconditioner.hpp"
#include "SingletonLogger.hpp"

namespace ProjectionSolvers {

	// helpers
	enum class StoppingCriterion { absolute, relative };

	using Preconditioner = Operator;

	inline void logInitialResidual(double r_0, Index max = 1) {
		auto& logger = SingletonLogger::instance();
		logger.buf << "start solving\n" 
		           << "|| r_" << std::setfill('0') << std::setw(std::to_string(max).length()) << '0' << " || = " << std::scientific << r_0;
		logger.log();
	}

	inline void logResidualReduction(double r_old, double r_new, Index i, Index max = 1) {
		auto& logger = SingletonLogger::instance();
		logger.buf << "|| r_" << std::setfill('0') << std::setw(std::to_string(max).length()) << i << " || = " << std::scientific << r_new << ", "
		           << "reduction factor = " << r_old / r_new;
		logger.log();
	}

	inline void logFinalResidual(double r_0, double r_n, Index i, Index max = 1) {
		auto& logger = SingletonLogger::instance();
		logger.buf << "stop solving\n"
		           << "|| r_" << std::setfill('0') << std::setw(std::to_string(max).length()) << i << " || = " << std::scientific << r_n << ", "
		           << "|| r_" << i << " || / || r_0 || = " << r_n / r_0;
		logger.log();
	}

	namespace Smoothers {

		inline std::vector<double> relaxedJacobi(
			AbstractPreconditioner<double>& A, // system matrix
			std::vector<double> const & b, // rhs
			boost::optional<std::vector<double>> const & x_0 = boost::none, // initial guess
			double omega = 1., // relaxation parameter
			Index n = 10000, // numb of iterations
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0 // log residual reduction on every i_log iteration (0 for never)
		) {
			// W (x_i+1 - x_i) = r_i,
			// A = L + D + U
			// W := omega^-1 diag(A)
			auto& logger = SingletonLogger::instance();
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
			     r_0 = b - A * x, r = r_0;
			decltype(x) r_new;
			double norm_r_0 = norm(r_0);
			logInitialResidual(norm_r_0, n);
			Index i;
			for (i = 1; i <= n; ++i) {
				x += omega * A.diagSubst(r);
				r_new = b - A * x;
				if (i_log && i % i_log == 0) logResidualReduction(norm(r), norm(r_new), i, n);
				r = r_new;
				if      (stop == StoppingCriterion::absolute && norm(r)            < eps) break;
				else if (stop == StoppingCriterion::relative && norm(r) / norm_r_0 < eps) break;
			}
			logFinalResidual(norm_r_0, norm(r), i > n ? n : i, n);
			if (eps != 0. && i > n) logger.wrn("Jacobi exceeded max numb of iterations");
			return x;
		}

		inline std::vector<double> forwSOR(
			AbstractPreconditioner<double>& A, // system matrix
			std::vector<double> const & b, // rhs
			boost::optional<std::vector<double>> const & x_0 = boost::none, // initial guess
			double omega = 1., // relaxation parameter
			Index n = 10000, // numb of iterations
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0 // log residual reduction on every i_log iteration (0 for never)
		) {
			// W (x_i+1 - x_i) = r_i,
			// A = L + D + R,
			// W := omega^-1 (omega L + D)
			// omega in (0, 2) for A = A^T > 0
			auto& logger = SingletonLogger::instance();
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
				r_0 = b - A * x, r = r_0;
			decltype(x) r_new;
			double norm_r_0 = norm(r_0);
			logInitialResidual(norm_r_0, n);
			Index i;
			for (i = 1; i <= n; ++i) {
				x += A.forwSubst(r, 1. / omega);
				r_new = b - A * x;
				if (i_log && i % i_log == 0) logResidualReduction(norm(r), norm(r_new), i, n);
				r = r_new;
				if      (stop == StoppingCriterion::absolute && norm(r)            < eps) break;
				else if (stop == StoppingCriterion::relative && norm(r) / norm_r_0 < eps) break;
			}
			logFinalResidual(norm_r_0, norm(r), i > n ? n : i, n);
			if (eps != 0. && i > n) logger.wrn("SOR exceeded max numb of iterations");
			return x;
		}

		inline std::vector<double> backSOR(
			AbstractPreconditioner<double>& A, // system matrix
			std::vector<double> const & b, // rhs
			boost::optional<std::vector<double>> const & x_0 = boost::none, // initial guess
			double omega = 1., // relaxation parameter
			Index n = 10000, // numb of iterations
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0 // log residual reduction on every i_log iteration (0 for never)
		) {
			// W (x_i+1 - x_i) = r_i,
			// A = L + D + R,
			// W := omega^-1 (omega U + D)
			// omega in (0, 2) for A = A^T > 0
			auto& logger = SingletonLogger::instance();
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
				r_0 = b - A * x, r = r_0;
			decltype(x) r_new;
			double norm_r_0 = norm(r_0);
			logInitialResidual(norm_r_0, n);
			Index i;
			for (i = 1; i <= n; ++i) {
				x += A.backSubst(r, 1. / omega);
				r_new = b - A * x;
				if (i_log && i % i_log == 0) logResidualReduction(norm(r), norm(r_new), i, n);
				r = r_new;
				if      (stop == StoppingCriterion::absolute && norm(r)            < eps) break;
				else if (stop == StoppingCriterion::relative && norm(r) / norm_r_0 < eps) break;
			}
			logFinalResidual(norm_r_0, norm(r), i > n ? n : i, n);
			if (eps != 0. && i > n) logger.wrn("SOR exceeded max numb of iterations");
			return x;
		}

		inline std::vector<double> SSOR(
			AbstractPreconditioner<double>& A, // system matrix
			std::vector<double> const & b, // rhs
			boost::optional<std::vector<double>> const & x_0 = boost::none, // initial guess
			double omega = 1., // relaxation parameter
			Index n = 10000, // numb of iterations
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0 // log residual reduction on every i_log iteration (0 for never)
		) {
			// W (x_i+1 - x_i) = r_i,
			// A = L + D + R,
			// W_0 := (2 - omega)^-1 (L + omega^-1 D) (omega^-1 D)^-1 (omega^-1 D + U)
			// omega in (0, 2) for A = A^T > 0
			auto& logger = SingletonLogger::instance();
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
				r_0 = b - A * x, r = r_0;
			decltype(x) r_new;
			double norm_r_0 = norm(r_0), omegaInv = 1. / omega;
			logInitialResidual(norm_r_0, n);
			Index i;
			for (i = 1; i <= n; ++i) {
				x += (2. * omegaInv - 1.) * A.backSubst(A.multDiag(A.forwSubst(r, omegaInv)), omegaInv);
				r_new = b - A * x;
				if (i_log && i % i_log == 0) logResidualReduction(norm(r), norm(r_new), i, n);
				r = r_new;
				if      (stop == StoppingCriterion::absolute && norm(r)            < eps) break;
				else if (stop == StoppingCriterion::relative && norm(r) / norm_r_0 < eps) break;
			}
			logFinalResidual(norm_r_0, norm(r), i > n ? n : i, n);
			if (eps != 0. && i > n) logger.wrn("SSOR exceeded max numb of iterations");
			return x;
		}

	}

	namespace Krylov {

		inline std::vector<double> CG(
			AbstractMultipliableMatrix<double> & A,
			std::vector<double> const & b,
			boost::optional<std::vector<double>> const & x_0 = boost::none,
			Index n = 0, // max numb of iters (0 for auto)
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0, // log residual reduction on every i_log iteration (0 for never)
			Index i_rec = 15 // recompute residual via b - A.x every i_rec iterations (0 for never)
		) {
			auto& logger = SingletonLogger::instance();
			// ini			
			if (!n) n = 3 * A.getOrder(); // max numb of iters
			Index i; // current iter
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
			     r = b - A * x,
			     p = r;
			auto r_x_r_0 = r * r, r_x_r = r_x_r_0;
			// dummy
			decltype(x) A_x_p;
			decltype(r_x_r) r_x_r_new, alpha;
			// start
			logInitialResidual(sqrt(r_x_r_0), n);
			for (i = 1; i <= n; ++i) {
				A_x_p = A * p;
				alpha = r_x_r / (A_x_p * p);
				x += alpha * p;
				r = !i_rec || i % i_rec ? r - alpha * A_x_p : b - A * x;
				r_x_r_new = r * r;
				if (i_log && i % i_log == 0) logResidualReduction(sqrt(r_x_r), sqrt(r_x_r_new), i, n);
				p = r + (r_x_r_new / r_x_r) * p;
				r_x_r = r_x_r_new;				
				if (stop == StoppingCriterion::absolute && sqrt(r_x_r)           < eps ||
				    stop == StoppingCriterion::relative && sqrt(r_x_r / r_x_r_0) < eps) break;
			}
			logFinalResidual(sqrt(r_x_r_0), sqrt(r_x_r), i > n ? n : i, n);
			if (i > n) logger.wrn("CG exceeded max numb of iterations");
			return x;
		}

		inline std::vector<double> PCG(
			Preconditioner const & B,
			AbstractMultipliableMatrix<double> & A,
			std::vector<double> const & b,
			boost::optional<std::vector<double>> const & x_0 = boost::none,
			Index n = 0,
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0,
			Index i_rec = 15
		) {
			auto& logger = SingletonLogger::instance();
			// ini
			if (!n) n = 3 * A.getOrder(); // max numb of iters
			Index i; // current iter
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
			     r = b - A * x,
			     z = B(r),
			     p = z;
			auto r_x_z_0 = r * z, r_x_z = r_x_z_0;
			// dummy
			decltype(x) A_x_p;
			double r_x_z_new, alpha;
			// start
			logger.log("here || . || denotes B-norm\n(B is a preconditioner and mimics inverse of A)");
			logInitialResidual(sqrt(r_x_z_0), n);
			for (i = 1; i <= n; ++i) {
				A_x_p = A * p;
				alpha = r_x_z / (A_x_p * p);
				x += alpha * p;
				r = !i_rec || i % i_rec ? r - alpha * A_x_p : b - A * x;
				z = B(r);
				r_x_z_new = r * z;
				if (i_log && i % i_log == 0) logResidualReduction(sqrt(r_x_z), sqrt(r_x_z_new), i, n);
				p = z + (r_x_z_new / r_x_z) * p;
				r_x_z = r_x_z_new;				
				if (stop == StoppingCriterion::absolute && sqrt(r_x_z)           < eps ||
				    stop == StoppingCriterion::relative && sqrt(r_x_z / r_x_z_0) < eps) break;
			}
			logFinalResidual(sqrt(r_x_z_0), sqrt(r_x_z), i > n ? n : i, n);
			if (i > n) logger.wrn("PCG exceeded max numb of iterations");
			return x;
		}

		inline std::vector<double> BiCGStab(
			AbstractMultipliableMatrix<double> & A,
			std::vector<double> const & b,
			boost::optional<std::vector<double>> const & x_0 = boost::none,
			Index n = 0,
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0,
			Index i_rec = 15
		) {
			auto& logger = SingletonLogger::instance();
			// ini
			if (!n) n = 3 * A.getOrder(); // max numb of iters
			Index i; // current iter
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
			     r = b - A * x, 
			     rbar = r, 
			     p = r;
			auto r_x_rbar = r * rbar,
			     norm_r_0 = norm(r), norm_r = norm_r_0;
			// dummy
			double norm_r_new, r_x_rbar_new, alpha, omega, beta;
			decltype(x) A_x_p, A_x_s, s;
			// start
			logInitialResidual(norm_r_0, n);
			for (i = 1; i <= n; ++i) {
				// (1) BiCG… step:
				A_x_p = A * p;
				alpha = r_x_rbar / (A_x_p * rbar);
				s = r - alpha * A_x_p;
				// s := residual after BiCG… step
				// check norm(s), if small enough: x += alpha * p and stop (Henk van der Vorst, p. 152)
				norm_r_new = norm(s);
				if (stop == StoppingCriterion::absolute && norm_r_new            < eps ||
					stop == StoppingCriterion::relative && norm_r_new / norm_r_0 < eps) {
					if (i_log && i % i_log == 0) logResidualReduction(norm_r, norm_r_new, i, n);
					logger.log("stopped after BiCG... iteration");
					norm_r = norm_r_new;
					x += alpha * p;
					break;
				}
				// (2) otherwise, make …Stab step:
				A_x_s = A * s;
				omega = (A_x_s * s) / (A_x_s * A_x_s);
				x += alpha * p + omega * s;
				r = i_rec && i % i_rec == 0 ? b - A * x : s - omega * A_x_s;
				r_x_rbar_new = r * rbar; 
				norm_r_new = norm(r);
				if (i_log && i % i_log == 0) logResidualReduction(norm_r, norm_r_new, i, n);
				beta = (r_x_rbar_new / r_x_rbar) * (alpha / omega);
				p = r + beta * (p - omega * A_x_p);
				r_x_rbar = r_x_rbar_new; 
				norm_r = norm_r_new;
				if (stop == StoppingCriterion::absolute && norm_r < eps ||
					stop == StoppingCriterion::relative && norm_r / norm_r_0 < eps) {
					logger.log("stopped after ...Stab iteration");
					break;
				}
			}
			logFinalResidual(norm_r_0, norm_r, i > n ? n : i, n);
			if (i > n) logger.wrn("BiCGStab exceeded max numb of iterations");
			return x;
		}

		inline std::vector<double> PBiCGStab(
			Preconditioner const & B,
			AbstractMultipliableMatrix<double> & A,
			std::vector<double> const & b,
			boost::optional<std::vector<double>> const & x_0 = boost::none,
			Index n = 0,
			double eps = 1e-12,
			StoppingCriterion stop = StoppingCriterion::absolute,
			Index i_log = 0,
			Index i_rec = 15
		) {
			auto& logger = SingletonLogger::instance();
			// ini
			if (!n) n = 3 * A.getOrder(); // max numb of iters
			Index i; // current iter
			auto x = x_0.value_or(std::vector<double>(A.getOrder())),
			     r = b - A * x, 
			     rbar = r, 
			     p = r;
			auto r_x_rbar = r * rbar,
			     norm_r_0 = norm(r), norm_r = norm_r_0;
			// dummy
			double norm_r_new, r_x_rbar_new, alpha, omega, beta;
			decltype(x) B_x_p, AB_x_p, B_x_s, AB_x_s, s;
			// start
			logInitialResidual(norm_r_0, n);
			for (i = 1; i <= n; ++i) {
				// (1) BiCG… step:
				B_x_p = B(p);
				AB_x_p = A * B_x_p;
				alpha = r_x_rbar / (AB_x_p * rbar);
				s = r - alpha * AB_x_p;
				// s := residual after BiCG… step
				// check norm(s), if small enough: x += alpha * p and stop (Henk van der Vorst, p. 152)
				norm_r_new = norm(s);
				if (stop == StoppingCriterion::absolute && norm_r_new            < eps ||
					stop == StoppingCriterion::relative && norm_r_new / norm_r_0 < eps) {
					if (i_log && i % i_log == 0) logResidualReduction(norm_r, norm_r_new, i, n);
					logger.log("stopped after BiCG... iteration");
					norm_r = norm_r_new;
					x += alpha * B_x_p;
					break;
				}
				// (2) otherwise, make …Stab step:
				B_x_s = B(s);
				AB_x_s = A * B_x_s;
				omega = (AB_x_s * s) / (AB_x_s * AB_x_s);
				x += alpha * B_x_p + omega * B_x_s;
				r = !i_rec || i % i_rec ? s - omega * AB_x_s : b - A * x;
				r_x_rbar_new = r * rbar; 
				norm_r_new = norm(r);
				if (i_log && i % i_log == 0) logResidualReduction(norm_r, norm_r_new, i, n);
				beta = (r_x_rbar_new / r_x_rbar) * (alpha / omega);
				p = r + beta * (p - omega * AB_x_p);
				r_x_rbar = r_x_rbar_new; 
				norm_r = norm_r_new;
				if (stop == StoppingCriterion::absolute && norm_r < eps ||
					stop == StoppingCriterion::relative && norm_r / norm_r_0 < eps) {
					logger.log("stopped after ...Stab iteration");
					break;
				}
			}
			logFinalResidual(norm_r_0, norm_r, i > n ? n : i, n);
			if (i > n) logger.wrn("BiCGStab exceeded max numb of iterations");
			return x;
		}

	}

}

//template <typename T>
//std::pair<std::vector<double>, size_t> BCG(T const &A,
//	std::vector<double> const &b,
//	std::vector<double> const &x0, 
//	double eps) {
//	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
//	std::vector<double> r1 = b - A*x0;
//	std::vector<double> r2 = r1;
//	std::vector<double> p1 = r1;
//	std::vector<double> p2 = r2;
//	std::vector<double> x = x0;
//	size_t iterCount = 1;
//	for (;iterCount <= 30000;iterCount++) {
//		double r1Xr2 = r1*r2;
//		std::vector<double> AXp1 = A*p1;
//		double a = r1Xr2 / (AXp1*p2);
//		x =x+ a*p1;
//		r1 = r1 - a*AXp1;
//		r2 = r2 - a*(A&p2);
//		double bet = (r1*r2) / r1Xr2;
//		if (abs(norm(r1) / norm(b)) < eps)
//			return std::pair<std::vector<double>, size_t>(x, iterCount);
//		if (bet == 0)
//			return std::pair<std::vector<double>, size_t>(x, -1);
//		p1 = r1 + bet*p1;
//		p2 = r2 + bet*p2;
//	}
//	return std::pair<std::vector<double>, size_t>(x, iterCount);
//}
//
//template <typename T>
//std::pair<std::vector<double>, size_t> ILUBCG(T &A,
//	std::vector<double> const &b,
//	std::vector<double> const &x0,
//	double eps) {
//	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
//	CSlRMatrix LU(A.ILU());
//	std::vector<double> r1 = b - A*x0;
//	std::vector<double> r2 = r1;
//	std::vector<double> p1 = r1;
//	std::vector<double> p2 = r2;
//	std::vector<double> x = x0;
//
//	size_t iterCount = 1;
//	for (;iterCount <= 30000;iterCount++) {
//		double r1Xr2 = r1*r2;
//		std::vector<double> AXp1 = LU.forwardSubstitution(A*(LU.backwardSubstitution(p1)));
//		double a = r1Xr2 / (AXp1*p2);
//		x = x + a*p1;
//		r1 = r1 - a*AXp1;
//		r2 = r2 - a*LU.forwardSubstitution(A&LU.backwardSubstitution(p2));
//		double bet = (r1*r2) / r1Xr2;
//		if (abs(norm(r1) / norm(b)) < eps)
//			return std::pair<std::vector<double>, size_t>(x, iterCount);
//		if (bet == 0)
//			return std::pair<std::vector<double>, size_t>(x, -1);
//		p1 = r1 + bet*p1;
//		p2 = r2 + bet*p2;
//	}
//	return std::pair<std::vector<double>, size_t>(x, iterCount);
//}
//
//template <typename T>
//std::pair<std::vector<double>, size_t> BiCGStab(T const &A,
//	std::vector<double> const &b,
//	std::vector<double> const &x0,
//	double eps) {
//	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
//	std::vector<double> r1 = b - A*x0;
//	std::vector<double> r2 = r1;
//	std::vector<double> p = r1;
//	std::vector<double> x = x0;
//	size_t iterCount = 1;
//	for (;iterCount <= 30000;iterCount++) {
//		double r1Xr2 = r1*r2;
//		std::vector<double> AXp = A*p;
//		double alpha = (r1Xr2) / (AXp*r2);
//		std::vector<double> s = r1 - alpha*AXp;
//		std::vector<double> AXs = A*s;
//		double omega = (AXs*s) / (AXs*AXs);
//		x = x + alpha*p + omega*s;
//		r1 = s - omega*AXs;
//		double beta = (r1*r2) / r1Xr2 * (alpha / omega);
//		if (abs(norm(r1) / norm(b)) < eps)
//			return std::pair<std::vector<double>, size_t>(x, iterCount);
//		if (beta == 0)
//			return std::pair<std::vector<double>, size_t>(x, -1);
//		p = r1 + beta*(p - omega*AXp);
//	}
//	return std::pair<std::vector<double>, size_t>(x, iterCount);
//}
//
//template <typename T>
//std::pair<std::vector<double>, size_t> ILUBiCGStab(T &A,
//	std::vector<double> const &b,
//	std::vector<double> const &x0,
//	double eps) {
//	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
//	CSlRMatrix LU(A.ILU());
//	std::vector<double> r1 = b - A*x0;
//	std::vector<double> r2 = r1;
//	std::vector<double> p = r1;
//	std::vector<double> x = x0;
//	size_t iterCount = 1;
//	for (;iterCount <= 30000;iterCount++) {
//		double r1Xr2 = r1*r2;
//		std::vector<double> AXp = LU.forwardSubstitution(A*(LU.backwardSubstitution(p)));
//		double alpha = (r1Xr2) / (AXp*r2);
//		std::vector<double> s = r1 - alpha*AXp;
//		std::vector<double> AXs = A*s;
//		double omega = (AXs*s) / (AXs*AXs);
//		x = x + alpha*p + omega*s;
//		r1 = s - omega*AXs;
//		double beta = (r1*r2) / r1Xr2 * (alpha / omega);
//		if (abs(norm(r1) / norm(b)) < eps)
//			return std::pair<std::vector<double>, size_t>(x, iterCount);
//		if (beta == 0)
//			return std::pair<std::vector<double>, size_t>(x, -1);
//		p = r1 + beta*(p - omega*AXp);
//	}
//	return std::pair<std::vector<double>, size_t>(x, iterCount);
//}