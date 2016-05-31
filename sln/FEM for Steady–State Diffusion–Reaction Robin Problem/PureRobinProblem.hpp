#pragma once
#include "unitCircle.hpp" // our mesh

// for any pure Robin problem:
inline double g_D(Node& p) {
	return 0.;
}
inline double g_N(Node& p) { // for any hom. problem (it is the case here)
	return 0.;
}

// (1)
inline double u(Node& p) { // a > 0
	return p.x() * p.x() + p.y() * p.y() + 1.;
}
inline double a(Node& p) { // a > 0
	return 7.;
}
inline double c(Node& p) { // c > c_0 >= 0
	return 2.;
}
inline double f(Node& p) {
	return 2 * p.x() * p.x() + 2 * p.y() * p.y() - 26.;
}
inline double kappa(Node& p) { // kappa > 0
	return 14. * (  1. / ( 1. + p.x() * p.x() + p.y() * p.y() ) - 1. );
}