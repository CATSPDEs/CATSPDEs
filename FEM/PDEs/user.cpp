#include "LinearEllipticPDE.hpp"

REAL force(Point const & point) {
	return 0.;
}

int main() {
	SymmetricMatrix I(3);
	I.identify();
	LinearEllipticPDE LaplacesEquation2D(I, force);
	std::cout << "div(D grad u) = f," << '\n'
		<< "D:" << '\n';
	LaplacesEquation2D.getDiffusionTensor().print();
	std::cout << '\n' << "f(0, 0):" << '\n'
			  << LaplacesEquation2D.getForceTerm(Point()) << std::endl;
	return 0;
}