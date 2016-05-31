#pragma once

class CurvilinearEdge {
	double _thetaStart, _thetaEnd;
	size_t _curveIndex;
public:
	CurvilinearEdge(double thetaStart, double thetaEnd, size_t curveIndex) 
		: _thetaStart(thetaStart)
		, _thetaEnd(thetaEnd)
		, _curveIndex(curveIndex) {
	}
	double thetaMiddle() { return (_thetaStart + _thetaEnd) / 2; }
	double& thetaStart() { return _thetaStart; }
	double& thetaEnd() { return _thetaEnd; }
	size_t& curveIndex() { return _curveIndex; }
};