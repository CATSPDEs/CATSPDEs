#pragma once
#include <stdexcept>
#include "Node.hpp"

class Triangle {
	Node* _p1; // counterclockwise! 
	Node* _p2;
	Node* _p3;
public:
	Triangle(Node* p1, Node* p2, Node* p3) : _p1(p1), _p2(p2), _p3(p3) {
		if (p1 == nullptr || p2 == nullptr || p3 == nullptr) throw std::logic_error("points must exist");
	}
	double area() { // norm of vector product underhood
		Node u = *_p2 - *_p1;
		Node v = *_p3 - *_p1;
		return u.crossProductNorm(v) / 2;
	}
};