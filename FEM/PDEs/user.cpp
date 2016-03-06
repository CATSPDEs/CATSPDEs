#include "LinearEllipticPDE.hpp"
#include "NonLinearEllipticPDE.hpp"

double force(Point const & point) {
	return 0.;
}

double f1(Point const & point) {
	return point.x() * point.x();
}

double f2(Point const & point) {
	return point.y() * point.y();
}

int main() {
	// laplace
	SymmetricMatrix I(3);
	LinearEllipticPDE LaplacesEquation3D(I.identify(), 0., force);
	std::cout << "div(D grad u) = f," << '\n'
		<< "D:" << '\n';
	LaplacesEquation3D.getDiffusionTensor().print();
	std::cout << "f(0, 0, 0):" << '\n'
			  << LaplacesEquation3D.getForceTerm(Point()) << std::endl;
	// my problem example
	SymmetricMatrixOfFunctions D(3);
	D(0, 0) = f1;
	D(1, 1) = f2;
	D(0, 1) = D(0, 2) = D(1, 2) = force;
	NonLinearEllipticPDE MyEqn2D(D, f1, force);
	std::cout << "Coeff of u_{xx} and u_{yy} at (1, 2):" << '\n'
			  << MyEqn2D.getDiffusionTensor(0, 0, Point(1., 2.)) << " " << MyEqn2D.getDiffusionTensor(1, 1, Point(1., 2.)) << '\n'
			  << "gamma(1, 2):" << '\n'
			  << MyEqn2D.getGamma(Point(1., 2.)) << std::endl;
	return 0;
}