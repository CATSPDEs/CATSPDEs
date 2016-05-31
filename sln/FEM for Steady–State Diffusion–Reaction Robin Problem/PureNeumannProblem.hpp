#pragma once
#include "unitCircle.hpp" // our mesh

// for any pure Neumann problem:
inline double g_D(Node& p) {
	return 0.;
}
inline double kappa(Node& p) { // kappa > 0
	return 0.;
}

// (1)
// u, c in P_1, a, f in P_2
//inline double u(Node& p) { // a > 0
//	return p.x() + 2. * p.y() + 3.;
//}
//inline double a(Node& p) { // a > 0
//	return p.x() * p.x() + p.x() * p.y();
//}
//inline double c(Node& p) { // c > c_0 >= 0
//	return 3. * p.y();
//}
//inline double f(Node& p) {
//	return p.x() * (3. * p.y() - 4.) + 2. * p.y() * (4. + 3. * p.y());
//}
//inline double g_N(Node& p) {
//	return p.x() * (p.x() + p.y()) * (p.x() + 2. * p.y());
//}

// (2)
// u, f in P_1, a in P_2, c = const
//inline double u(Node& p) { // a > 0
//	return p.x() + 2. * p.y() + 3.;
//}
//inline double a(Node& p) { // a > 0
//	return p.x() * p.x() + p.x() * p.y();
//}
//inline double c(Node& p) { // c > c_0 >= 0
//	return 1.;
//}
//inline double f(Node& p) {
//	return p.y() - 3. * ( p.x() - 1. );
//}
//inline double g_N(Node& p) {
//	return p.x() * ( p.x() + p.y() ) * ( p.x() + 2. * p.y() );
//}

// (3)
// u, a, f in P_1, c = const 
//inline double u(Node& p) { // a > 0
//	return p.x() + 2. * p.y() + 3.;
//}
//inline double a(Node& p) { // a > 0
//	return p.x() + p.y();
//}
//inline double c(Node& p) { // c > c_0 >= 0
//	return 1.;
//}
//inline double f(Node& p) {
//	return p.x() + 2. * p.y();
//}
//inline double g_N(Node& p) {
//	return (p.x() + p.y()) * (p.x() + 2. * p.y());
//}

// (4)
// u in P_2, a = const, c = const AND g_N = const on bndry
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
inline double g_N(Node& p) {
	return 14. * ( p.x() * p.x() + p.y() * p.y() );
}

// (5)
// u in P_2, a = const, c = const BUT g_N != const on bndry
//inline double u(Node& p) { // a > 0
//	return (p.x() - .5) * (p.x() - .5) + 4. * p.y() * p.y() + 1.;
//}
//inline double a(Node& p) { // a > 0
//	return 7.;
//}
//inline double c(Node& p) { // c > c_0 >= 0
//	return 2.;
//}
//inline double f(Node& p) {
//	return 2* ( (p.x() - .5) * (p.x() - .5) + 4. * p.y() * p.y() - 34. );
//}
//inline double g_N(Node& p) {
//	return 7. * p.x() * ( 2. * p.x() - 1. ) + 56. * p.y() * p.y();
//}