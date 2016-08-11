#pragma once
#include "IMultipliable.hpp"

template <typename T>
class ITransposeMultipliable : public IMultipliable<T> {
private:
	bool _t;
	virtual vector<T> _multByTranspose(vector<T> const &) const = 0;
public:
	ITransposeMultipliable() : _t(false) {}
	ITransposeMultipliable& t() {
		_t = !_t;
		return *this;
	}
	vector<T> operator*(vector<T> const & b) {
		bool dummy = _t;
		_t = false;
		return dummy ? _multByTranspose(b) : _mult(b);
	}
};