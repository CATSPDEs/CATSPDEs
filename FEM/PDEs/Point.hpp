#pragma once
#include "real.hpp"

class Point {
	REAL _x, _y, _z;
public:
	explicit Point(REAL x = 0., REAL y = 0., REAL z = 0.) : _x(x), _y(y), _z(z) {}
	Point& add(Point const & point) {
		_x += point._x;
		_y += point._y;
		_z += point._z;
		return *this;
	}
	Point& mult(REAL scaler) {
		_x *= scaler;
		_y *= scaler;
		_z *= scaler;
		return *this;
	}
};