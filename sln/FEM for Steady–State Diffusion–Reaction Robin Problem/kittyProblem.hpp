#pragma once
#include "kitty.hpp" // our mesh

// eqn

inline double a(Node& p) { // a > 0
	return 1.;
}
inline double c(Node& p) {
	return 0.;
}
inline double f(Node& p) {
	return 0.;
}

// BCs

inline double g_D(Node& p) {
	return p.x() * p.x() + p.y() * p.y();
}
DirichletBC Dirichlet(g_D); // pure Neumann everywhere
vector<AbstractBC*> ListOfBCs = { &Dirichlet };