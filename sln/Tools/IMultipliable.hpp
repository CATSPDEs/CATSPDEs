#pragma once
#include <vector>

using namespace std;

template <typename T>
class IMultipliable {
protected:
	virtual vector<T> _mult(vector<T> const &) = 0;
public:
	virtual ~IMultipliable() {}
	virtual vector<T> operator*(vector<T> const & b) {
		return _mult(b);
	}
};