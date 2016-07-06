#pragma once
#include<vector>
#include<type_traits>
#include"IRealMatrix.hpp"

template <typename T>
std::vector<double> CG(T const &A, 
	std::vector<double> const &b, 
	std::vector<double> const &x0, 
	double const eps) {
	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
	std::vector<double> r = b - A*x0;
	std::vector<double> p = r;
	std::vector<double> x = x0;
	for (unsigned maxIter = 0; /*maxIter < 10000 &&*/ abs(norm(r) / norm(b)) > eps; ++maxIter) {
		double  rXr = r*r;
		std::vector<double> AXp = (A*p);
		double a = rXr / (AXp*p);
		x = x + a*p;
		r = r - a*AXp;
		double b = (r*r) / rXr;
		p = r + b*p;
	}
	return x;
}

template <typename T>
std::pair<std::vector<double>, size_t> BCG(T const &A,
	std::vector<double> const &b,
	std::vector<double> const &x0, 
	double const eps) {
	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
	std::vector<double> r1 = b - A*x0;
	std::vector<double> r2 = r1;
	std::vector<double> p1 = r1;
	std::vector<double> p2 = r2;
	std::vector<double> x = x0;
	size_t iterCount = 1;
	for (;iterCount <= 30000;iterCount++) {
		double r1Xr2 = r1*r2;
		std::vector<double> AXp1 = A*p1;
		double a = r1Xr2 / (AXp1*p2);
		x =x+ a*p1;
		r1 = r1 - a*AXp1;
		r2 = r2 - a*(A&p2);
		double bet = (r1*r2) / r1Xr2;
		if (abs(norm(r1) / norm(b)) < eps)
			return std::pair<std::vector<double>, size_t>(x, iterCount);
		if (bet == 0)
			return std::pair<std::vector<double>, size_t>(x, -1);
		p1 = r1 + bet*p1;
		p2 = r2 + bet*p2;
	}
	return std::pair<std::vector<double>, size_t>(x, iterCount);
}

template <typename T>
std::pair<std::vector<double>, size_t> ILUBCG(T &A,
	std::vector<double> const &b,
	std::vector<double> const &x0,
	double const eps) {
	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
	CSlRMatrix LU(A.ILU());
	std::vector<double> r1 = b - A*x0;
	std::vector<double> r2 = r1;
	std::vector<double> p1 = r1;
	std::vector<double> p2 = r2;
	std::vector<double> x = x0;

	size_t iterCount = 1;
	for (;iterCount <= 30000;iterCount++) {
		double r1Xr2 = r1*r2;
		std::vector<double> AXp1 = LU.forwardSubstitution(A*(LU.backwardSubstitution(p1)));
		double a = r1Xr2 / (AXp1*p2);
		x = x + a*p1;
		r1 = r1 - a*AXp1;
		r2 = r2 - a*LU.forwardSubstitution(A&LU.backwardSubstitution(p2));
		double bet = (r1*r2) / r1Xr2;
		if (abs(norm(r1) / norm(b)) < eps)
			return std::pair<std::vector<double>, size_t>(x, iterCount);
		if (bet == 0)
			return std::pair<std::vector<double>, size_t>(x, -1);
		p1 = r1 + bet*p1;
		p2 = r2 + bet*p2;
	}
	return std::pair<std::vector<double>, size_t>(x, iterCount);
}

template <typename T>
std::pair<std::vector<double>, size_t> BiCGStab(T const &A,
	std::vector<double> const &b,
	std::vector<double> const &x0,
	double const eps) {
	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
	std::vector<double> r1 = b - A*x0;
	std::vector<double> r2 = r1;
	std::vector<double> p = r1;
	std::vector<double> x = x0;
	size_t iterCount = 1;
	for (;iterCount <= 30000;iterCount++) {
		double r1Xr2 = r1*r2;
		std::vector<double> AXp = A*p;
		double alpha = (r1Xr2) / (AXp*r2);
		std::vector<double> s = r1 - alpha*AXp;
		std::vector<double> AXs = A*s;
		double omega = (AXs*s) / (AXs*AXs);
		x = x + alpha*p + omega*s;
		r1 = s - omega*AXs;
		double beta = (r1*r2) / r1Xr2 * (alpha / omega);
		if (abs(norm(r1) / norm(b)) < eps)
			return std::pair<std::vector<double>, size_t>(x, iterCount);
		if (beta == 0)
			return std::pair<std::vector<double>, size_t>(x, -1);
		p = r1 + beta*(p - omega*AXp);
	}
	return std::pair<std::vector<double>, size_t>(x, iterCount);
}

template <typename T>
std::pair<std::vector<double>, size_t> ILUBiCGStab(T &A,
	std::vector<double> const &b,
	std::vector<double> const &x0,
	double const eps) {
	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
	CSlRMatrix LU(A.ILU());
	std::vector<double> r1 = b - A*x0;
	std::vector<double> r2 = r1;
	std::vector<double> p = r1;
	std::vector<double> x = x0;
	size_t iterCount = 1;
	for (;iterCount <= 30000;iterCount++) {
		double r1Xr2 = r1*r2;
		std::vector<double> AXp = LU.forwardSubstitution(A*(LU.backwardSubstitution(p)));
		double alpha = (r1Xr2) / (AXp*r2);
		std::vector<double> s = r1 - alpha*AXp;
		std::vector<double> AXs = A*s;
		double omega = (AXs*s) / (AXs*AXs);
		x = x + alpha*p + omega*s;
		r1 = s - omega*AXs;
		double beta = (r1*r2) / r1Xr2 * (alpha / omega);
		if (abs(norm(r1) / norm(b)) < eps)
			return std::pair<std::vector<double>, size_t>(x, iterCount);
		if (beta == 0)
			return std::pair<std::vector<double>, size_t>(x, -1);
		p = r1 + beta*(p - omega*AXp);
	}
	return std::pair<std::vector<double>, size_t>(x, iterCount);
}