#pragma once
#include "Node.hpp"

typedef bool(*Predicate)(Node&); // function pointer

inline bool constTrue (Node&) { return true; }
inline bool constFalse(Node&) { return false; }