#pragma once
#include "AbstractQuadratureRule.hpp"

// quad rule for master node
class MasterNodeQuadratureRule
	: public AbstractQuadratureRule<0> {
	// singleton
	MasterNodeQuadratureRule() : AbstractQuadratureRule { 
		{ {} }, { 1. }
	}
	{}
	MasterNodeQuadratureRule(MasterNodeQuadratureRule const &);
	MasterNodeQuadratureRule& operator=(MasterNodeQuadratureRule const &);
public:
	static auto& instance() {
		static MasterNodeQuadratureRule single;
		return single;
	}
};