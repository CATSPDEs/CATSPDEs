#pragma once
#include <functional>
#include "Node.hpp"

template <LocalIndex N, LocalIndex M>
using Mapping = std::function<Node<M>(Node<N> const &)>; // general mapping : R^N —> R^M

template <LocalIndex N, LocalIndex M>
inline auto constZero(Node<N> const &) { 
	Node<M> res;
	if (M == 1) res = 0.;
	else res.fill(0.);
	return res; 
}

template <LocalIndex N, LocalIndex M>
inline auto constUnity(Node<N> const &) {
	Node<M> res;
	if (M == 1) res = 1.;
	else res.fill(1.);
	return res;
}

template <LocalIndex N>
using VectorField = Mapping<N, N>;

template <LocalIndex N>
using ScalarField = Mapping<N, 1>;
