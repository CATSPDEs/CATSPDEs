#pragma once
#include "Function.hpp"

class BC {
	Function _D, _N, _RU, _R;
	// u = D on Г1,
	// u_{n} = N on Г2,
	// u_{n} + R0 u = R on Г3 
public:
	explicit BC(Function D = nullptr, Function N = nullptr,
				Function R = nullptr, Function R0 = oneFunc)
				: _D(D), _N(N), _R(R), _RU(R0) {}
	double dirichlet(Point const &p) const { return _D(p); }
	double neumann(Point const &p) const { return _N(p); }
	double robinU(Point const &p) const { return _RU(p); }
	double robin(Point const &p) const { return _R(p); }
};
