#pragma once

// exact soln

inline double u(Node& p, double t) {
	return p.x() * p.x() + p.y() * p.y() + 1.;
}

// IBCs

inline double kappa(Node& p, double t) {
	return 0.;
}

inline double g_N(Node& p, double t) {
	return 14. * (p.x() * p.x() + p.y() * p.y());
}

inline double g_D(Node& p, double t) {
	return 0.;
}

inline double initialVelocity(Node& p) {
	return 14. * (p.x() * p.x() + p.y() * p.y());
}

// PDE

inline double chi(Node& p) {
	return 0.;
}

inline double sigma(Node& p) { 
	return p.x() * p.x() + p.y() * p.y() + 1.;
}

inline double a(Node& p) {
	return 7.;
}

inline double f(Node& p, double t) {
	return 2 * p.x() * p.x() + 2 * p.y() * p.y() - 26.;
}