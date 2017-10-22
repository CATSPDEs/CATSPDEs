#pragma once
#include <fstream>
#include <set>
#include "Element.hpp"
#include "Parameters.hpp"
#include "vector.hpp"
#include "Predicate.hpp"

/*
	Alexander Žilyakov, Oct 2016
*/

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element

template <LocalIndex D, LocalIndex N>
class AbstractMesh {
protected:
	std::vector<Node<D>>              _nodes; // nodes (i.e. P-matrix),
	std::vector<std::array<Index, N>> _elements;
	// vector of arrays of indicies of nodes that span an element of the meash at hand 
	// e.g triangle, rectangle, simplex etc.
	// (i.e. T-matrix or connectivity matrix)
	std::set<Index> _ghostElements;
public:
	AbstractMesh() {}
	AbstractMesh(std::vector<Node<D>> const & nodes, std::vector<std::array<Index, N>> const & elements)
		: _nodes(nodes), _elements(elements)
	{}
	virtual ~AbstractMesh() {}
	// pure virtual member funcs
	virtual AbstractMesh& refine(Index numbOfRefinements = 1) = 0;
	virtual AbstractMesh& import(std::istream& from = std::cin) = 0;
	virtual void          export(std::ostream& to = std::cout, Parameters const & param = {}) const = 0;
	// get numb of nodes
	Index numbOfNodes() const { return _nodes.size(); }
	// get numb of elements
	virtual Index numbOfElements() const { return _elements.size(); }
	// get ith node
	Node<D> getNode(Index i) const { return _nodes[i]; }
	auto const & getNodes() const { return _nodes; }
	// get nodes’ indicies of eth element
	virtual std::array<Index, N> getNodesIndicies(Index e) const {
		return _elements[e];
	}
	// get eth element
	virtual Element<D, N> getElement(Index e) const {
		Element<D, N> el;
		for (LocalIndex i = 0; i < N; ++i)
			el[i] = _nodes[_elements[e][i]];
		return el;
	}
	// get indicies of new elements (added after refinement) 
	// contained in COARSE element e or
	virtual std::vector<Index> getFineNeighborsIndicies(Index e) const = 0;
	// i/o
	AbstractMesh& import(std::string const & fromStr) {
		std::ifstream from(fromStr);
		return import(from);
	}
	void export(std::string const & toStr, Parameters const & param = {}) const {
		std::ofstream to(toStr);
		export(to, param);
	}
};