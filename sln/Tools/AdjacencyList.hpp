#pragma once
#include <vector>
#include <set>

// useful data structure which is used to 
// generate pattern of sparse (CSC-like) matrix
// from the FEM mesh
typedef std::vector<std::set<size_t>> AdjacencyList;