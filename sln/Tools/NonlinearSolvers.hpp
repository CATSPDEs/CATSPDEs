#pragma once
#include <boost/optional/optional.hpp>
#include <boost/circular_buffer.hpp>
#include "DenseMatrix.hpp"
#include "Mapping.hpp"

namespace NonlinearSolvers {

	struct FixedPointData {
		Operator f; // residual operator
		Operator g; // iteration operator
	};

	inline void logFixedPoint(double r, Index i, Index max = 1) {
		auto& logger = SingletonLogger::instance();
		auto& w = std::setw(std::to_string(max).length());
		logger.buf << "|| f(x_" << std::setfill('0') << w << i << ") || = " << std::scientific << r;
		logger.log();
	}

	inline std::vector<double> FixedPointMethod(
		FixedPointData const & fp,
		std::vector<double> const & x_0,
		Index n = 100, // max numb of iters
		double eps = 1e-7
	) {
		auto& logger = SingletonLogger::instance();
		logger.log("start Fixed Point method");
		auto x = x_0;
		auto r_norm = norm(fp.f(x));
		logFixedPoint(r_norm, 0, n);
		if (r_norm < eps) {
			logger.log("stop Fixed Point method");
			return x;
		}
		decltype(x) x_trial, g;
		decltype(r_norm) r_norm_new, omega;
		Index i;
		for (i = 1; i <= n; ++i) {
			omega = 1.;
			r_norm_new = r_norm + 1.;
			g = fp.g(x);
			while (r_norm_new > r_norm && omega > 1e-10) {
				x_trial = omega * g + (1. - omega) * x;
				r_norm_new = norm(fp.f(x_trial));
				omega /= 2.;
			}
			x = x_trial;
			r_norm = r_norm_new;
			logFixedPoint(r_norm, i, n);
			if (r_norm < eps) break;
		}
		logger.log("stop Fixed Point method");
		if (i > n) logger.wrn("Fixed Point Method exceeded max numb of iterations");
		return x;
	}
	
	inline std::vector<double> AndersonMixingMethod(
		FixedPointData const & fp,
		std::vector<double> const & x_0,
		Index m = 0,
		Index n = 100, // max numb of iters
		double eps = 1e-7
		) {
		auto& logger = SingletonLogger::instance();
		if (m == 0) {
			logger.log("m = 0 -> Fixed Point Method");
			return FixedPointMethod(fp, x_0, n, eps);
		}
		boost::circular_buffer<std::vector<double>> F(m), G(m + 1);
		logger.log("start Fixed Point method");
		// i = 0
		auto x = x_0;
		auto r = fp.f(x);
		auto r_norm = norm(r);
		logFixedPoint(r_norm, 0, n);
		if (r_norm < eps) {
			logger.log("stop Fixed Point method");
			return x;
		}
		G.push_front(fp.g(x));
		// i = 1
		x = G[0];
		auto r_new = fp.f(x);
		r_norm = norm(r_new);
		logFixedPoint(r_norm, 1, n);
		if (r_norm < eps) {
			logger.log("stop Fixed Point method");
			return x;
		}
		G.push_front(fp.g(x));
		F.push_front(r_new - r);
		r = r_new;
		logger.log("stop Fixed Point method");
		// i = 2, 3, . . .
		logger.log("start Anderson Mixing method");
		Index i;
		for (i = 2; i <= n; ++i) {
			m = F.size();
			DenseMatrix<double> A(m);
			std::vector<double> b(m);
			for (Index k = 0; k < m; ++k) {
				b[k] = F[k] * r;
				for (Index j = 0; j < m; ++j)
					A(k, j) = F[k] * F[j];
			}
			auto beta = A.GaussElimination(b);
			decltype(beta) alpha(m + 1);
			alpha[0] = 1. - beta[0];
			for (Index k = 1; k < m; ++k) alpha[k] = beta[k - 1] - beta[k];
			alpha[m] = beta[m - 1];
			std::fill(x.begin(), x.end(), 0.);
			for (Index k = 0; k < m + 1; ++k) x += alpha[k] * G[k];
			r_new = fp.f(x);
			r_norm = norm(r_new);
			logFixedPoint(r_norm, i, n);
			if (r_norm < eps) break;
			G.push_front(fp.g(x));
			F.push_front(r_new - r);
			r = r_new;
		}
		logger.log("stop Anderson Mixing method");
		if (i > n) logger.wrn("Anderson Mixing Method exceeded max numb of iterations");
		return x;
	}

}