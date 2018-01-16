#pragma once
#include "Mapping.hpp"
#include "Predicate.hpp"

template <LocalIndex N, LocalIndex M>
class BoundaryCondition {
	std::vector<Mapping<N, M>>  _f;
	Predicate<N>                _p;
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

	template <LocalIndex M>
	using BoundaryCondition2D = BoundaryCondition<2, M>;

		using ScalarBoundaryCondition2D = BoundaryCondition2D<1>;
		using VectorBoundaryCondition2D = BoundaryCondition2D<2>;