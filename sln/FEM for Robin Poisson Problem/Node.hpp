#pragma once
#include <cmath>
#include <iostream>

using namespace std;

class Node {
	double _x;
	double _y;
public:
	explicit Node(double x = 0., double y = 0.) : _x(x), _y(y) {}
	double& x() { return _x; }
	double& y() { return _y; }
	Node& operator+=(Node const & p) {
		_x += p._x;
		_y += p._y;
		return *this;
	}
	Node operator+(Node const & p) const {
		return Node(*this) += p;
	}
	Node& operator-=(Node const & p) { 
		_x -= p._x;
		_y -= p._y;
		return *this;
	}
	Node operator-(Node const & p) const {
		return Node(*this) -= p;
	}
	double operator*(Node const & p) { // dot product
		return _x * p._x + _y * p._y;
	}
	double crossProductNorm(Node const & p) {
		return _x * p._y - _y * p._x;
	}
	double norm() {
		return sqrt((*this) * (*this));
	}
	friend ostream& operator<<(ostream& out, Node const & p) {
		return out << p._x << ' ' << p._y << '\n';
	}
	friend istream& operator>>(istream& in, Node& p) {
		return in >> p._x >> p._y;
	}
};