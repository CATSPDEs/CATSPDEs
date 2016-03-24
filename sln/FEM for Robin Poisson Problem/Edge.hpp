#pragma once
#include <stdexcept>
#include "Node.hpp"

class Edge {
	Node* _beg;
	Node* _end;
public:
	Edge(Node* beg, Node* end) : _beg(beg), _end(end) {
		if (beg == nullptr || end == nullptr) throw std::logic_error("points must exist");
	}
	double length() { // compute length of Edge
		return (*_beg - *_end).norm();
	}
};
