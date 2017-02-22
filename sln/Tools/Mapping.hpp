#pragma once
#include <functional>
#include "Node.hpp"

// general mapping f : R^N —> R^M
template <LocalIndex N, LocalIndex M>
using Mapping = std::function<Node<M>(Node<N> const &)>;

	template <LocalIndex N>
	using VectorField = Mapping<N, N>;

	using VectorField2D = VectorField<2>;
	using VectorField3D = VectorField<3>;

	template <LocalIndex N>
	using ScalarField = Mapping<N, 1>;

	using ScalarField2D = ScalarField<2>;
	using ScalarField3D = ScalarField<3>;

	template <LocalIndex M>
	using Curve = Mapping<1, M>;

	using Curve2D = Curve<2>;
	using Curve3D = Curve<3>;
