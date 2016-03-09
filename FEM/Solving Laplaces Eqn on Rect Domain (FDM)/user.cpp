#include "LinearEllipticPDE.hpp"
#include "BC.hpp"

double g(Point const & point) { // u = g on Г
	return point.x() * point.y();
}

int main() {
	// Laplace's eqn on the plane
	ConstTensor I(2, 0.);
	I(0, 0) = I(1, 1) = 1.;
	LinearEllipticPDE LaplacesEquation2D(I, 0., zeroFunc);
	BC myBC(g); // simple Dirichlet problem
	std::cout << "u(2, 3) = " << myBC.dirichlet(Point(2., 3.)) << std::endl;
	return 0;
}