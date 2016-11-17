#pragma once
#include <functional>
#include "Node.hpp"

// D := node dimension

template <LocalIndex D>
using Predicate = std::function<bool(Node<D> const &)>;

template <LocalIndex D>
inline bool constTrue (Node<D> const &) { return true; }

template <LocalIndex D>
inline bool constFalse(Node<D> const &) { return false; }