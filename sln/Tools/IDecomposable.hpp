#pragma once
#include <vector>

template <typename T>
class IDecomposable { // interface for 
public:
	virtual ~IDecomposable() {}
	virtual std::vector<T> forwardSubstitution()  const = 0;
	virtual std::vector<T> backwardSubstitution() const = 0;
};