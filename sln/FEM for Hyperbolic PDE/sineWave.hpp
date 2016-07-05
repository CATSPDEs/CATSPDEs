#pragma once

// model

inline double u(Node& p, double t) {
	return cos(5. * (p.x() + p.y() + t));
}

// eqn

inline double chi(Node& p) {
	return 1.;
}
inline double sigma(Node& p) {
	return .1;
}
inline double a(Node& p) {
	return p.x() * p.y();
}
inline double f(Node& p, double t) {
	return 25. * (2. * p.x() * p.y() - 1.) * cos(5. * (t + p.x() + p.y())) +
	       5. * (p.x() + p.y() - .1) * sin(5. * (t + p.x() + p.y()));
}

// ICs

inline double initialPosition(Node& p) {
	return u(p, 0.);
}
inline double initialVelocity(Node& p) {
	return -5. * sin( 5. * ( p.x() + p.y() ) );
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
	return -5. * p.x() * p.y() * sin(5. * (t + p.x() + p.y()));;
}

// (3) Robin boundary
inline double G3_R(Node& p, double t) {
	return 1.;
}
inline double G3_N(Node& p, double t) {
	return cos(5. * (t + p.x() + p.y())) + 5. * p.x() * p.y() * sin(5. * (t + p.x() + p.y()));
}
