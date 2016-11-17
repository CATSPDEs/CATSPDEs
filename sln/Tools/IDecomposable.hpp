#pragma once
#include <vector>

template <typename T>
class IDecomposable { // interface for 
public:
	virtual ~IDecomposable() {}
	virtual IDecomposable& decompose() = 0; // compute (incomplete) LDU (LU, LDL^T etc.) decomposition
	virtual vector<T> forwardSubstitution()  const = 0; // LU.x = b, (1) solve L.y = b
	virtual vector<T> backwardSubstitution() const = 0; //           (2) solve U.x = y
};