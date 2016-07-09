#pragma once
#include "IMultipliable.hpp"

template <typename T>
class ITransposeMultipliable : public IMultipliable<T> {
private:
	bool _t;
	virtual std::vector<T> _multByTranspose(std::vector<T> const &) const = 0;
public:
	ITransposeMultipliable() : _t(false) {}
	ITransposeMultipliable& t() {
		_t = !_t;
		return *this;
	}
	std::vector<T> operator*(std::vector<T> const & b) {
		auto v = _t ? _multByTranspose(b) : _mult(b);
		_t = false;
		return v;
	}
};