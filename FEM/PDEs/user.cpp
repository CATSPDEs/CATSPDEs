#include "LinearEllipticPDE.hpp"

REAL force(Point const & point) {
	return 0.;
}

int main() {
	SymmetricMatrix I(2);
	try {
		std::cout << I.identify();
	}
	catch (std::out_of_range& e) {
		std::cout << e.what();
	}
	//LinearEllipticPDE LaplacesEquation2D(I, force);
	//Point dummy(0.);
	//std::cout << "div(D grad u) = f," << '\n'
	//	<< "D:" << '\n'
	//	<< LaplacesEquation2D.getDiffusionTensor()
	//	<< "f(0, 0):" << '\n';
		//<< LaplacesEquation2D.getForceTerm(dummy);
	return 0;
}