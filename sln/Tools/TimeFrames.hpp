#pragma once
#include "vector.hpp"

class TimeFrames {
private:
	vector<double> _t;
public:
	explicit TimeFrames(double initialTime = 0., double endTime = 1., unsigned refCount = 0) 
		: _t(pow(2, refCount) + 1)
		{
		if (initialTime < 0.) throw invalid_argument("t0 < 0 — time cannot take negative values!");
		if (initialTime >= endTime) throw invalid_argument("we should have t0 < t1 < … < tn");
		size_t stepsCount = _t.size() - 1;
		double step = (endTime - initialTime) / stepsCount;
		for (size_t m = 0; m <= stepsCount; ++m)
			_t[m] = initialTime + m * step;
	}
	TimeFrames& refine(unsigned refCount = 1) {
		double initialTime = _t[0], endTime = _t[_t.size() - 1];
		_t.resize((_t.size() - 1) * pow(2, refCount) + 1);
		size_t stepsCount = _t.size() - 1;
		double step = (endTime - initialTime) / stepsCount;
		for (size_t m = 0; m <= stepsCount; ++m)
			_t[m] = initialTime + m * step;
		return *this;
	}
	double at(size_t i) const {
		return _t.at(i);
	}
	double operator[](size_t i) const {
		return _t[i];
	}
	size_t size() const {
		return _t.size();
	}
	ostream& save(ostream& output = cout) {
		return output << _t;
	}
};