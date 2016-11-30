#pragma once
#include <functional>
#include "Node.hpp"

template <LocalIndex N, LocalIndex M>
using Mapping = std::function<Node<M>(Node<N> const &)>; // general mapping : R^N —> R^M

	template <LocalIndex N>
	using VectorField = Mapping<N, N>;

		using VectorField2D = VectorField<2>;

	template <LocalIndex N>
	using ScalarField = Mapping<N, 1>;

		using ScalarField2D = ScalarField<2>;
