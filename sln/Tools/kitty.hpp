#pragma once
#include <fstream>
#include "Triangulation.hpp"

ifstream kittyNodes("../Tools/kittyNodes.dat"),
         kittyTriangles("../Tools/kittyTriangles.dat");
Triangulation Omega(kittyNodes, kittyTriangles);
