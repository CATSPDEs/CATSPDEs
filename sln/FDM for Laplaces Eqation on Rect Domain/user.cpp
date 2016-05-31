#include "SecondOrderPDE.hpp"
#include "BC.hpp"

double gamma(Point const & p) {
	return p.x() * p.x();
}

double D(Point const & p) { // u = g on Г
	return p.x() * p.y();
}

int main() {
	SecondOpderPDE HelmholtzEqn(oneFunc, gamma, zeroFunc);
	BC myBC(D); // simple Dirichlet problem
	std::cout << "u(2, 3) = " << myBC.dirichlet(Point(2., 3.)) << std::endl;
	return 0;
}