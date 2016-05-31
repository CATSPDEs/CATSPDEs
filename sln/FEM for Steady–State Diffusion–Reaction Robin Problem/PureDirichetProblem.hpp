#pragma once
#include "unitCircle.hpp" // our mesh
//#include "unitSquare.hpp"

// for any pure Dirichlet problem:
inline double g_N(Node& p) {
	return 0.;
}
inline double g_D(Node& p) {
	return u(p);
}
inline double kappa(Node& p) { // kappa > 0
	return 10e+50;
}

// (1)
// u, c in P_1, a, f in P_2
inline double u(Node& p) { // a > 0
	return p.x() + 2. * p.y() + 3.;
}
inline double a(Node& p) { // a > 0
	return p.x() * p.x() + p.x() * p.y();
}
inline double c(Node& p) { // c > c_0 >= 0
	return 3. * p.y();
}
inline double f(Node& p) {
	return p.x() * (3. * p.y() - 4.) + 2. * p.y() * (4. + 3. * p.y());
}