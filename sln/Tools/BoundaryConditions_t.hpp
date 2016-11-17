#pragma once
#include <stdexcept>
#include <vector>
#include <memory>
#include "Function_t.hpp"
#include "Predicate.hpp"

class AbstractBC_t { // boundary condition (any kind)
protected:
	Predicate _validAt; // if _validAt(node), then at node this BCs are applied
	// –n[r] . (a[r] nabla u[r, t]) = _R[r, t] (u[r, t] – _D[r, t]) – _N[r, t]
public:
	explicit AbstractBC_t(Predicate validAt = constTrue)
		: _validAt(validAt) {}
	virtual ~AbstractBC_t() {}
	bool validAt(Node& p) const { return _validAt(p); }
	// we will redefine these later for Dirichlet, Neumann, and Robin BCs
	virtual double DirichletCondition(Node& p, double t) const = 0;
	virtual double NeumannValue      (Node& p, double t) const = 0;
	virtual double RobinCoefficient  (Node& p, double t) const = 0;
};

// more specific:

class DirichletBC_t : public AbstractBC_t {
	Function_t _D;
	// u[r, t] – _D[r, t] = –n[r] . (a[r] nabla u[r, t]) / 10e50 ~ 0 or
	// u[r, t] ~ _D[r, t]
	// we approximate Dirichlet BCs w/ Robin BCs here
public:
	explicit DirichletBC_t(Function_t D = zeroFunc_t, Predicate validAt = constTrue)
		: AbstractBC_t(validAt)
		, _D(D) {}
	// hom. BCs by default
	// u[r, t] ~ 0
	double DirichletCondition(Node& p, double t) const { return _D(p, t); }
	double NeumannValue      (Node& p, double t) const { return 0.; }
	double RobinCoefficient  (Node& p, double t) const { return 10e50; }
};

class NeumannBC_t : public AbstractBC_t {
	Function_t _N;
	// n[r] . (a[r] nabla u[r, t]) = _N[r, t]
public:
	explicit NeumannBC_t(Function_t N = zeroFunc_t, Predicate validAt = constTrue)
		: AbstractBC_t(validAt)
		, _N(N) {}
	// hom. BCs by default
	// n[r] . (a[r] nabla u[r, t]) = 0
	double DirichletCondition(Node& p, double t) const { return 0.; }
	double NeumannValue      (Node& p, double t) const { return _N(p, t); }
	double RobinCoefficient  (Node& p, double t) const { return 0.; }
};

class RobinBC_t : public AbstractBC_t {
	Function_t _R, _N;
	// n[r] . (a[r] nabla u[r, t]) + _R[r, t] u[r, t] = _N[r, t]
public:
	explicit RobinBC_t(Function_t R = oneFunc_t, Function_t N = zeroFunc_t, Predicate validAt = constTrue)
		: AbstractBC_t(validAt)
		, _R(R)
		, _N(N) {}
	// hom. BCs by default
	// n[r] . (a[r] nabla u[r, t]) + u[r, t] = 0
	double DirichletCondition(Node& p, double t) const { return 0.; }
	double NeumannValue      (Node& p, double t) const { return _N(p, t); }
	double RobinCoefficient  (Node& p, double t) const { return _R(p, t); }
};

// list of BCs (any kind)
class BoundaryConditions_t {
	vector<shared_ptr<AbstractBC_t>> _BCs_t;
	size_t _current; // index of current BCs to apply
public:
	BoundaryConditions_t(vector<shared_ptr<AbstractBC_t>> const & BCs_t) : _BCs_t(BCs_t), _current(0) {
		if (BCs_t.size() == 0) throw invalid_argument("list of BCs cannot be empty");
	}
	void defineBCsAt(Node& p) {
		for (size_t i = 0; i < _BCs_t.size(); ++i)
			if (_BCs_t[i]->validAt(p)) {
				_current = i;
				return;
			}
	}
	double DirichletCondition(Node& p, double t) const { return _BCs_t[_current]->DirichletCondition(p, t); }
	double NeumannValue      (Node& p, double t) const { return _BCs_t[_current]->NeumannValue(p, t); }
	double RobinCoefficient  (Node& p, double t) const { return _BCs_t[_current]->RobinCoefficient(p, t); }
};
