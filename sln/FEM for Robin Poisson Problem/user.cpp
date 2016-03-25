#include <iostream>
#include "Triangulation.hpp"

int main() {
	Triangulation K(Node(0, 0), Node(4, 2), .9);
	cout << K.area(0) << ' ' << K.length(0, 0) << '\n'; 
	// area of the very 1st triangle and length of its 1st side
	// (i.e. side against its 1st node)
	// it happens to be hypotenuse (and legs measure is equal to unity)
	return 0;
}