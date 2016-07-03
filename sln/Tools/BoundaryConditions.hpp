#pragma once
#include <stdexcept>
#include <vector>
#include "Function.hpp"
#include "Predicate.hpp"

using namespace std;

class AbstractBC { // boundary condition (any kind)
protected:
	Predicate _validAt; // if _validAt(node), then at node this BCs are applied
						// –n . (a nabla u) = _R (u – _D) – _N
public:
	explicit AbstractBC(Predicate validAt = constTrue)
		: _validAt(validAt) {}
	virtual ~AbstractBC() {}
	bool validAt(Node& p) const { return _validAt(p); }
	// we will redefine these later for Dirichlet, Neumann, and Robin BCs
	virtual double DirichletCondition(Node& p) const = 0;
	virtual double NeumannValue      (Node& p) const = 0;
	virtual double RobinCoefficient  (Node& p) const = 0;
};

// more specific:

class DirichletBC : public AbstractBC {
	Function _D;
	// u – _D = –n . (a nabla u) / 10e50 ~ 0 or
	// u ~ _D
	// we approximate Dirichlet BCs w/ Robin BCs here
public:
	explicit DirichletBC(Function D = zeroFunc, Predicate validAt = constTrue)
		: AbstractBC(validAt)
		, _D(D) {}
	// hom. BCs by default
	// u ~ 0
	double DirichletCondition(Node& p) const { return _D(p); }
	double NeumannValue      (Node& p) const { return 0.; }
	double RobinCoefficient  (Node& p) const { return 10e50; }
};

class NeumannBC : public AbstractBC {
	Function _N;
	// n . (a nabla u) = _N
public:
	explicit NeumannBC(Function N = zeroFunc, Predicate validAt = constTrue)
		: AbstractBC(validAt)
		, _N(N) {}
	// hom. BCs by default
	// n . (a nabla u) = 0
	double DirichletCondition(Node& p) const { return 0.; }
	double NeumannValue      (Node& p) const { return _N(p); }
	double RobinCoefficient  (Node& p) const { return 0.; }
};

class RobinBC : public AbstractBC {
	Function _R, _N;
	// n . (a nabla u) + _R u = _N
public:
	explicit RobinBC(Function R = oneFunc, Function N = zeroFunc, Predicate validAt = constTrue)
		: AbstractBC(validAt)
		, _R(R)
		, _N(N) {}
	// hom. BCs by default
	// n . (a nabla u) + u = 0
	double DirichletCondition(Node& p) const { return 0.; }
	double NeumannValue      (Node& p) const { return _N(p); }
	double RobinCoefficient  (Node& p) const { return _R(p); }
};

// list of BCs (any kind)
class BoundaryConditions {
	vector<AbstractBC*> _BCs;
	size_t _current; // index of current BCs to apply
public:
	BoundaryConditions(vector<AbstractBC*> const & BCs) : _BCs(BCs), _current(0) {
		if (BCs.size() == 0) throw invalid_argument("list of BCs cannot be empty");
	}
	void defineBCsAt(Node& p) {
		for (size_t i = 0; i < _BCs.size(); ++i)
			if (_BCs[i]->validAt(p)) {
				_current = i;
				return;
			}
	}
	double DirichletCondition(Node& p) const { return _BCs[_current]->DirichletCondition(p); }
	double NeumannValue      (Node& p) const { return _BCs[_current]->NeumannValue(p); }
	double RobinCoefficient  (Node& p) const { return _BCs[_current]->RobinCoefficient(p); }
};