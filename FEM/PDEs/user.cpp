#include "LinearEllipticPDE.hpp"

double f1(Point const & point) {
	return point.x() * point.x();
}

double f2(Point const & point) {
	return point.y() * point.y();
}

int main() {
	// laplace
	ConstTensor I(3, 0.);
	I(0, 0) = I(1, 1) = I(2, 2) = 1.;
	LinearEllipticPDE LaplacesEquation3D(I, 0., zeroFunc);
	std::cout << "div(D grad u) = f, D:" << '\n';
	LaplacesEquation3D.diffusionTensor().save(std::cout);
	std::cout << "f(0, 0, 0): " << LaplacesEquation3D.forceTerm(Point()) << std::endl;
	// my problem example
	Tensor D(2);
	//D.load(std::cin);
	D(0, 0) = f1;
	D(1, 1) = f2;
	D(0, 1) = zeroFunc;
	NonLinearEllipticPDE MyEqn2D(D, f1, zeroFunc);
	Point dummy(1., 2.);
	std::cout << "Coeff of u_{xx}, u_{yy} etc. at (1, 2):" << '\n';
	MyEqn2D.diffusionTensor(dummy).save(std::cout);
	std::cout << "gamma(1, 2): " << MyEqn2D.gamma(dummy) << '\n'
			  << "f(1, 2): " << MyEqn2D.forceTerm(dummy) << std::endl;
	return 0;
}