#pragma once
#include "Mapping.hpp"
#include "Predicate.hpp"

template <LocalIndex N, LocalIndex M>
class BoundaryCondition {
	std::vector<Mapping<N, M>>  _f;
	Predicate<N>                 _p;
public:
	explicit BoundaryCondition(Mapping<N, M> const & f, Predicate<N> const & p = constTrue<N>)
		: _f({ f })
		, _p(p) {}
	explicit BoundaryCondition(std::vector<Mapping<N, M>> const & f, Predicate<N> const & p = constTrue<N>)
		: _f(f)
		, _p(p) {}
	Node<M> operator()        (Node<N> const & p, Index i = 0) const { return (_f[i])(p); }
	bool    shouldBeEnforcedAt(Node<N> const & p) const { return _p(p); }
};

template <LocalIndex N>
using ScalarBoundaryCondition = BoundaryCondition<N, 1>;

using ScalarBoundaryCondition1D = ScalarBoundaryCondition<1>;
using ScalarBoundaryCondition2D = ScalarBoundaryCondition<2>;
using ScalarBoundaryCondition3D = ScalarBoundaryCondition<3>;

template <LocalIndex N>
using VectorBoundaryCondition = BoundaryCondition<N, N>;

using VectorBoundaryCondition1D = VectorBoundaryCondition<1>;
using VectorBoundaryCondition2D = VectorBoundaryCondition<2>;
using VectorBoundaryCondition3D = VectorBoundaryCondition<3>;
