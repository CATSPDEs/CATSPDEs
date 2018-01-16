#pragma once
#include <numeric> // inner_product
#include "array.hpp"

// D := dimension of the node (1 for segments, 2 for triangulations etc.)

template <LocalIndex D> 
using Node = SmartArray<double, D>;

	using Node2D = Node<2>;
	using Node3D = Node<3>;
