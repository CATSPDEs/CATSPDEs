#pragma once
#include <vector>
#include <set>

typedef std::vector<std::set<size_t>> AdjacencyList;

// useful data structure which is used to 
// generate portrait of sparse (CRS-like) matrix
// from the FEM mesh