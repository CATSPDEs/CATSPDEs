#pragma once
#include<vector>
#include<type_traits>
#include"IRealMatrix.hpp"

template <typename T>
std::vector<double> CG(T const &A, std::vector<double> const &b, std::vector<double> const &x0, double const eps) {
	static_assert(std::is_base_of<IRealMatrix, T>::value, "matrix class must be a descendant of IRealMatrix");
	std::vector<double> r = b - A*x0;
	std::vector<double>p = r;
	std::vector<double> x = x0;
	while (abs(norm(r) / norm(b)) > eps) {
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
