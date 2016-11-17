#pragma once
#include <boost/tuple/tuple.hpp>
#include <boost/optional/optional.hpp>
#include "AbstractMultipliableMatrix.hpp"

extern SingletonLogger& logger;

namespace ProjectionSolvers {

	namespace Smoothers {

		inline
		boost::tuple<
			std::vector<double>, // soln
			std::vector<double> // norms of residuals
		> 
		Jacobi(
			AbstractMultipliableMatrix<double>& A, // system matrix
			std::vector<double> const & b, // rhs
			boost::optional<std::vector<double>> const & x_0 = boost::none, // initial guess
			double omega = 1., // relaxation parameter
			Index numbOfIterations = 10000, // numb of iterations
			double const eps = 10e-17
		) {
			// W (x_i+1 - x_i) = r_i,
			// W := omega^-1 diag(A)
			auto x = x_0.value_or(std::vector<double>(A.getOrder(), 0.)),
			     r = b - A * x;
			std::vector<double> residualsNorms;
			for (Index i = 0; i < numbOfIterations; ++i) {
				residualsNorms.emplace_back(norm(r));
				for (Index j = 0; j < A.getOrder(); ++j)
					x[j] += omega * r[j] / A(j, j);
				r = b - A * x;
				if (norm(r) < eps) break;
			}
			return boost::make_tuple(x, residualsNorms);
		}

		inline
		boost::tuple<
			std::vector<double>, // soln
			std::vector<double> // norms of residuals
		> 
		GaussSeidel(
			SymmetricCSlCMatrix<double>& A, // system matrix
			std::vector<double> const & b, // rhs
			boost::optional<std::vector<double>> const & x_0 = boost::none, // initial guess
			double omega = 1., // relaxation parameter
			Index numbOfIterations = 10000, // numb of iterations
			double const eps = 10e-17
		) {
			// W (x_i+1 - x_i) = r_i,
			// A = L + D + R,
			// W := omega^-1 (omega L + D)
			// omega in (0, 2) for A = A^T > 0
			auto x = x_0.value_or(std::vector<double>(A.getOrder(), 0.)),
			     r = b - A * x;
			std::vector<double> residualsNorms;
			for (Index i = 0; i < numbOfIterations; ++i) {
				residualsNorms.emplace_back(norm(r));
				if (norm(r) < eps) break;
				auto y(r);
				for (Index j = 0; j < A.getOrder(); ++j) {
					y[j] *= omega / A(j, j);
					x[j] += y[j];
					for (Index k = A._colptr[j]; k < A._colptr[j + 1]; ++k)
						y[A._rowind[k]] -= A._lval[k] * y[j];
				}
				r = b - A * x;
			}
			return boost::make_tuple(x, residualsNorms);
		}

	}

	namespace Krylov {

		inline
		boost::tuple<
			std::vector<double>, // soln vector
			std::vector<double> // norms of residuals
			//Index // numb of iterations
		> CG(
			AbstractMultipliableMatrix<double> & A,
			std::vector<double> const & b,
			boost::optional<std::vector<double>> const & x_0 = boost::none,
			double const eps = 10e-17
		) {
			auto x = x_0.value_or(std::vector<double>(A.getOrder(), 0.)),
			     r = b - A * x,
			     p = r;
			auto r_x_r = r * r, normOfb = norm(b);
			Index i, maxNumbOfIterations = 1.5 * A.getOrder();
			std::vector<double> residualsNorms;
			for (i = 0; i < maxNumbOfIterations; ++i) {
				residualsNorms.emplace_back(sqrt(r_x_r));
				if (residualsNorms.back() < eps) break;
				if (i % 50 == 0) {
					logger.buf << "|| r_" << i << " || = " << std::scientific << residualsNorms.back();
					logger.log();
				}
				auto A_x_p = A * p;
				auto alpha = r_x_r / (A_x_p * p);
				x += alpha * p;
				r -= alpha * A_x_p;
				auto r_x_r_new = r * r;
				p = r + (r_x_r_new / r_x_r) * p;
				r_x_r = r_x_r_new;
			}
			if (i == maxNumbOfIterations)
				logger.wrn(
					"CG exceeded max numb of iterations: i = " + std::to_string(maxNumbOfIterations) + ", matrix size = " + std::to_string(A.getOrder()) + '\n' +
					"|| r_i || = " + std::to_string(norm(r))
				);
			return boost::make_tuple(x, residualsNorms);
		}

		inline
		boost::tuple<
			std::vector<double>, // soln vector
			std::vector<double> // norms of residuals
		> PCG(
			AbstractMultipliableMatrix<double> & A,
			std::vector<double> const & b,
			std::function<std::vector<double>(
				std::vector<double> const &
			)> const & B, // preconditioner 
			boost::optional<std::vector<double>> const & x_0 = boost::none,
			double const eps = 10e-17
		) {
			auto x = x_0.value_or(std::vector<double>(A.getOrder(), 0.)),
			     r = b - A * x,
			     z = B(r),
			     p = z;
			auto r_x_z = r * z;
			decltype(x) A_x_p;
			decltype(r_x_z) r_x_z_new, alpha;
			Index i, maxNumbOfIterations = 1.5 * A.getOrder();
			std::vector<double> residualsNorms;
			for (i = 0; i < maxNumbOfIterations; ++i) {
				residualsNorms.emplace_back(norm(r));
				logger.buf << "|| r_" << i << " || = " << std::scientific << residualsNorms.back();
				logger.log();
				if (residualsNorms.back() < eps) break;
				A_x_p = A * p;
				alpha = r_x_z / (A_x_p * p);
				x += alpha * p;
				r = b - A * x; // r -= alpha * A_x_p;
				z = B(r);
				r_x_z_new = r * z;
				p = z + (r_x_z_new / r_x_z) * p;
				r_x_z = r_x_z_new;
			}
			if (i == maxNumbOfIterations)
				logger.wrn(
					"CG exceeded max numb of iterations: i = " + std::to_string(maxNumbOfIterations) + ", matrix size = " + std::to_string(A.getOrder()) + '\n' +
					"|| r_i || = " + std::to_string(norm(r))
				);
			return boost::make_tuple(x, residualsNorms);
		}

	}

}

//template <typename T>
//std::pair<std::vector<double>, size_t> BCG(T const &A,
//	std::vector<double> const &b,
//	std::vector<double> const &x0, 
//	double const eps) {
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
//	double const eps) {
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
//	double const eps) {
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
//	double const eps) {
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