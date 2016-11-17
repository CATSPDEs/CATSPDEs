#pragma once
#include <functional>
#include "Node.hpp"

template <LocalIndex D>
using Function = std::function<double(Node<D> const &)>;

template <LocalIndex D>
inline double emptyFunc(Node<D> const &) { return 0.; } // for default parameters

template <LocalIndex D>
inline double zeroFunc(Node<D> const &) { return 0.; }

template <LocalIndex D>
inline double oneFunc(Node<D> const &) { return 1.; }
