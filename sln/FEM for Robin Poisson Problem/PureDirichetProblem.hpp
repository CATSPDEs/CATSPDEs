#pragma once
#include "Triangulation.hpp" // our mesh

inline double a(Node& p) { // a > 0
	return p.x() * p.x() + 7. * p.x() * p.y() + 2. * p.y() + 3.;
}

inline double c(Node& p) { // c > c_0 >= 0
	return 0.;
}

inline double f(Node& p) {
	return -25. * p.x() - 14. * p.y() - 6.;
}

inline double g_D(Node& p) {
	return 2. * p.x() + 3. * p.y() + 1.;
}

inline double g_N(Node& p) {
	return 0.;
}

inline double kappa(Node& p) { // kappa > 0
	return 10e+50;
}

Triangulation Omega(Node(-1., -1.), Node(1., 1.), .5); // simple square mesh