#pragma once
#include <boost/optional/optional.hpp>
#include "Mapping.hpp"

namespace NonlinearSolvers {

	inline void logPicard(double r, double d, Index i, Index max = 1) {
		auto& logger = SingletonLogger::instance();
		auto& w = std::setw(std::to_string(max).length());
		logger.buf 
			<< "|| r_" << std::setfill('0') << w << i - 1 << " || = " << std::scientific << r << ",\n"
			<< "|| x_" << std::setfill('0') << w << i << " - x_" << std::setfill('0') << w << i - 1 << " || = " << d;
		logger.log();
	}

	using PicardOperator = std::function<std::vector<double>(std::vector<double> const &, double&)>;

	inline std::vector<double> PicardIteration(
		PicardOperator const & phi,
		std::vector<double> const & x_0,
		double omega = 1.,
		Index n = 100, // max numb of iters
		double eps = 1e-7
	) {
		// x_k+1 = phi(x_k) until convergence
		auto x_new = x_0, x_old = x_0;
		double delta, r_norm;
		auto& logger = SingletonLogger::instance();
		logger.log("start Picard iteration");
		Index i;
		for (i = 1; i <= n; ++i) {
			x_new = omega * phi(x_old, r_norm) + (1. - omega) * x_old;
			delta = norm(x_new - x_old);
			logPicard(r_norm, delta, i, n);
			if (r_norm < eps) break;
			x_old = x_new;
		}
		logger.log("stop Picard iteration");
		if (i > n) logger.wrn("Picard exceeded max numb of iterations");
		return x_new;
	}

}