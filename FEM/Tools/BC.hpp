#pragma once
#include "Function.hpp"

class BC {
	Function _D, _N, _R0, _R1, _R;
	// u = D on Г1,
	// u_{n} = N on Г2,
	// R0 u + R1 u_{n} = R on Г3 
public:
	explicit BC(Function D = nullptr, Function N = nullptr,
				Function R = nullptr, Function R0 = oneFunc, Function R1 = oneFunc)
				: _D(D), _N(N), _R(R), _R0(R0), _R1(R1) {}
	double dirichlet(Point const &p) const { return _D(p); }
	double neumann(Point const &p) const { return _N(p); }
	double robinU(Point const &p) const { return _R0(p); }
	double robinUn(Point const &p) const { return _R1(p); }
	double robin(Point const &p) const { return _R(p); }
};
