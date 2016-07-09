#pragma once
#include <vector>

template <typename T>
class IMultipliable {
protected:
	virtual std::vector<T> _mult(std::vector<T> const &) const = 0;
public:
	virtual ~IMultipliable() {}
	virtual std::vector<T> operator*(std::vector<T> const & b) {
		return _mult(b);
	}
};