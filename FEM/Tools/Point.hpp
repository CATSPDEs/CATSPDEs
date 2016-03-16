#pragma once

class Point {
	double _x, _y, _z;
public:
	explicit Point(double x = 0., double y = 0., double z = 0.) : _x(x), _y(y), _z(z) {}
	Point& add(Point const & point) {
		_x += point._x;
		_y += point._y;
		_z += point._z;
		return *this;
	}
	Point& mult(double scaler) {
		_x *= scaler;
		_y *= scaler;
		_z *= scaler;
		return *this;
	}
	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
};

/*
class Point2D : public Point {
protected:
	double _y;
public:
	explicit Point2D(double x = 0., double y = 0.) : Point(x), _y(y) {}
	double y() const { return _y; }
};
*/

/*
class Point {
	double _x, _y, _z;
public:
	explicit Point(double x = 0., double y = 0., double z = 0.) : _x(x), _y(y), _z(z) {}
	Point& add(Point const & point) {
		_x += point._x;
		_y += point._y;
		_z += point._z;
		return *this;
	}
	Point& mult(double scaler) {
		_x *= scaler;
		_y *= scaler;
		_z *= scaler;
		return *this;
	}
	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
};
*/