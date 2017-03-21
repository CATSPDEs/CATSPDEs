#pragma once
#include "AbstractFiniteElement.hpp"

/*
	Alexander Žilyakov, March 2017
*/

template <LocalIndex D, LocalIndex N, LocalIndex M>
class FEInterpolant {
	std::vector<double> _DOFs; // degrees of freedom
	// aggregation
	AbstractFiniteElement<D, N, M> const * _FE;
	AbstractMesh<D, N> const * _mesh;
public:
	FEInterpolant(std::vector<double> const & DOFs, AbstractFiniteElement<D, N, M> const & FE, AbstractMesh<D, N> const & mesh)
		: _DOFs(DOFs)
		, _FE(&FE)
		, _mesh(&mesh) {
		if (FE.numbOfDOFs(mesh) != DOFs.size()) throw std::invalid_argument("invalid numb of DOFs");
	}
	FEInterpolant(
		ScalarField<D> const & f, AbstractFiniteElement<D, N, M> const & FE, AbstractMesh<D, N> const & mesh, 
		boost::optional<Index&> activeElementIndex = boost::none
	)
		: _DOFs(FE.numbOfDOFs(mesh))
		, _FE(&FE)
		, _mesh(&mesh) {
		// index of the active element
		Index activeElementIndexScoped;
		Index& e = activeElementIndex.value_or(activeElementIndexScoped);
		for (e = 0; e < mesh.numbOfElements(); ++e) {
			auto nodes = FE.getDOFsNodes(mesh, e);
			auto indicies = FE.getDOFsNumeration(mesh, e);
			for (Index i = 0; i < indicies.size(); ++i)
				_DOFs[indicies[i]] = f(nodes[i]);
		}
	}
	auto& DOFs() { return _DOFs; }
	Node<M> operator()(Node<D> const & p, SignedIndex elementIndex = -1) const {
		if (elementIndex == -1) {
			if (D == 1) { // segment mesh
				auto const & nodes = _mesh->getNodes();
				if (nodes.back() < p || p < nodes.front()) throw std::invalid_argument("invalid domain");
				if (p == nodes.back()) elementIndex = _mesh->numbOfElements() - 1;
				else {
					// binary search for element
					auto upper = std::upper_bound(nodes.begin(), nodes.end(), p);
					elementIndex = upper - nodes.begin() - 1;
				}
			}
			else throw std::logic_error("not implemented");
		}
		// shapes
		auto shapes = _FE->getShapesOf(_mesh->getElement(elementIndex));
		Index n = shapes.size();
		std::vector<Node<M>> sValues(n);
		std::transform(shapes.begin(), shapes.end(), sValues.begin(), [&](auto const & s) {
			return s(p);
		});
		// weights
		auto indicies = _FE->getDOFsNumeration(*_mesh, elementIndex);
		std::vector<double> weights;
		weights.reserve(indicies.size());
		for (Index i : indicies)
			weights.emplace_back(_DOFs[i]);
		// result
		return weights * sValues;
	}
	Node<D * M> grad(Node<D> const & p, SignedIndex elementIndex = -1) const {
		if (elementIndex == -1) throw std::logic_error("not implemented");
		// shapes gradients
		auto sGrads = _FE->getSGradsOf(_mesh->getElement(elementIndex));
		Index n = sGrads.size();
		std::vector<Node<2 * M>> gValues(n);
		std::transform(sGrads.begin(), sGrads.end(), gValues.begin(), [&](auto const & g) {
			return g(p);
		});
		// weights
		auto indicies = _FE->getDOFsNumeration(*_mesh, elementIndex);
		std::vector<double> weights;
		weights.reserve(indicies.size());
		for (Index i : indicies)
			weights.emplace_back(_DOFs[i]);
		// result
		return weights * gValues;
	}
};
	
	using SegmentFEInterpolant = FEInterpolant<1, 2, 1>;
	
	using TriangularScalarFEInterpolant = FEInterpolant<2, 3, 1>;