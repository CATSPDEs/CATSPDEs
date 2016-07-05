#pragma once

// model

inline double u(Node& p, double t) {
	return p.x() - 2. * p.y() + t * t;
}

// eqn

inline double chi(Node& p) {
	return 2.;
}
inline double sigma(Node& p) {
	return 1.;
}
inline double a(Node& p) {
	return 1.;
}
inline double f(Node& p, double t) {
	return 2. * ( t + 2.);
}

// ICs

inline double initialPosition(Node& p) {
	return u(p, 0.);
}
inline double initialVelocity(Node& p) {
	return 0.;
}

// BCs

// (1) Dirichlet boundary
inline bool G1(Node& p) {
	if (p.y() == 0. || p.y() == 1.) return true;
	return false;
}
inline double G1_D(Node& p, double t) {
	return u(p, t); 
}

// (2) Neumann boundary
inline bool G2(Node& p) {
	if (p.x() == 1.) return true;
	return false;
}
inline double G2_N(Node& p, double t) {
	return 1.;
}

// (3) Robin boundary
inline bool G3(Node& p) {
	if (p.x() == 0.) return true;
	return false;
}
inline double G3_R(Node& p, double t) {
	return 1.;
}
inline double G3_N(Node& p, double t) {
	return -1. + t * t + p.x() - 2. * p.y();
}

