#pragma once
#include <unordered_map>
#include <string>
#include "Mapping.hpp"

template <LocalIndex N, LocalIndex M>
class SmartMapping {
	// functor for general SmartMapping f : R^N —> R^M
	// represented by std::function (e.g. by lambda)
	Mapping<N, M> _f;
	std::unordered_map<Node<N>, Node<M>> _savedImages;
	Index _capacity = 1000;
public:
	template <typename T, typename Requires = decltype(Mapping<N, M>(std::declval<T&&>()))>
	SmartMapping(T f) : _f(std::move(f)) {}
	Index capacity() const { return _capacity; }
	Index& capacity() { return _capacity; }
	Node<M> operator()(Node<N> const & p) {
		auto kvp = _savedImages.find(p);
		if (kvp != _savedImages.end()) return kvp->second;
		if (_savedImages.size() >= _capacity) throw std::runtime_error("numb of images needed to be saved > SmartMapping capacity = " + std::to_string(_capacity));
		auto val = _f(p);
		_savedImages[p] = val;
		return val;
	}
	// map _f on several nodes
	std::vector<Node<M>> operator()(std::vector<Node<N>> const & nodes) {
		std::vector<Node<M>> images(nodes.size());
		std::transform(nodes.begin(), nodes.end(), images.begin(), [this](auto const & p) {
			return operator()(p);
		});
		return images;
	}
	SmartMapping& saveImagesOf(std::vector<Node<N>> const & nodes) {
		for (auto const & p : nodes) _savedImages[p] = _f(p);
		if (_savedImages.size() > _capacity) _capacity = _savedImages.size();
		return *this;
	}
};

	template <LocalIndex N>
	using SmartVectorField = SmartMapping<N, N>;

	using SmartVectorField2D = SmartVectorField<2>;
	using SmartVectorField3D = SmartVectorField<3>;

	template <LocalIndex N>
	using SmartScalarField = SmartMapping<N, 1>;

	using SmartScalarField1D = SmartScalarField<1>;
	using SmartScalarField2D = SmartScalarField<2>;
	using SmartScalarField3D = SmartScalarField<3>;

template <LocalIndex N, LocalIndex M>
std::vector<SmartMapping<N, M>> smartify(std::vector<Mapping<N, M>> const & mappings) {
	std::vector<SmartMapping<N, M>> smartMappings;
	smartMappings.reserve(mappings.size());
	for (auto const & map : mappings) smartMappings.push_back(map);
	return smartMappings;
}