#pragma once
#include "Node.hpp"

template <LocalIndex N, LocalIndex M>
class AbstractMapping { 
	// abstract functor for general mapping f : R^N —> R^M
public:
	virtual Node<M> operator()(Node<N> const &) const = 0;
	virtual ~AbstractMapping() {}
};
