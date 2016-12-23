#pragma once
#include <unordered_map>
#include "AbstractMapping.hpp" // !!!
#include "Mapping.hpp"

template <LocalIndex N, LocalIndex M>
class SmartMapping : public AbstractMapping<N, M> {
	// functor for general SmartMapping f : R^N —> R^M
	// represented by std::function (e.g. by lambda)
	Mapping<N, M> _f;
	std::unordered_map<Node<N>, Node<M>> _savedImages;
public:
	//SmartMapping(
	//	Mapping<N, M> const & f,
	//	std::unordered_map<Node<N>, Node<M>> const & savedImages
	//) : _f(f), _savedImages(savedImages) {}
	template <typename T, typename Requires = decltype(Mapping<N, M>(std::declval<T&&>()))>
	SmartMapping(T f) : _f(std::move(f)) {}
	Node<M> operator()(Node<N> const & p) const final {
		auto kvp = _savedImages.find(p);
		return kvp != _savedImages.end() ? kvp->second : _f(p);
	}
	// map _f on several nodes
	std::vector<Node<M>> operator()(std::vector<Node<N>> const & nodes) const {
		std::vector<Node<M>> images(nodes.size());
		std::transform(nodes.begin(), nodes.end(), images.begin(), [this](auto const & p) {
			return operator()(p);
		});
		return images;
	}
	SmartMapping& saveImagesOf(std::vector<Node<N>> const & nodes) {
		for (auto const & p : nodes) _savedImages[p] = _f(p);
		return *this;
	}
	//// scalar / dot product multiplication of SmartMappings
	//SmartMapping<N, 1> operator*(SmartMapping const & m) {
	//	// (1) construct function
	//	auto f = [&](Node<N> const & p) {
	//		return _f(p) * m._f(p);
	//	};
	//	// (2) recompute saved images in a smart way
	//	std::unordered_map<Node<N>, double> savedImages;
	//	for (auto const & kvp :   _savedImages) savedImages[kvp.first] = kvp.second            * m(kvp.first);
	//	for (auto const & kvp : m._savedImages) savedImages[kvp.first] = operator()(kvp.first) * kvp.second;
	//	return { f, savedImages };
	//}
	//friend SmartMapping operator*(AbstractMultipliableMatrix<double>& A, SmartMapping const & m) {
	//	if (M != A.numbOfRows() || M != A.numbOfCols()) throw std::invalid_argument("only square operators are supported");
	//	// (1) construct function
	//	auto f = [&](Node<N> const & p) {
	//		return A * m._f(p);
	//	};
	//	// (2) recompute saved images in a smart way
	//	auto savedImages = m._savedImages;
	//	for (auto& kvp : savedImages) kvp.second = A * kvp.second;
	//	return { f, savedImages };
	//}
};

	template <LocalIndex N>
	using SmartVectorField = SmartMapping<N, N>;

	using SmartVectorField2D = SmartVectorField<2>;
	using SmartVectorField3D = SmartVectorField<3>;

	template <LocalIndex N>
	using SmartScalarField = SmartMapping<N, 1>;

	using SmartScalarField2D = SmartScalarField<2>;
	using SmartScalarField3D = SmartScalarField<3>;