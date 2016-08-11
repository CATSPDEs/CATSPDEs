#include "StokesMatrix.hpp"
#include <algorithm>

vector<double> StokesMatrix::_mult(vector<double> const & u){ // compute v = S.u
	vector<double> v(_w, 0.);
	size_t i, n = _A11.getOrder(),
	          m = _B1.numbOfRows();
	// A11 & B1
	vector<double> vdummy = _A11 * u;
	copy(vdummy.begin(), vdummy.end(), v.begin());
	vdummy = _B1 * u;
	copy(vdummy.begin(), vdummy.begin() + _B1.numbOfRows(), v.end() - _B1.numbOfRows);

	// B1^T
	//copy(u.end() - _B1.numbOfRows(), u.end(), dummy.begin());
	//dummy = B1.t() * dummy;
	//for (i = 0; i < A11.getOrder(); ++i)
	//	v[i] += dummy[i];

	// A22
	vector<double> udummy(_A11.getOrder, 0);
	copy(u.begin() + _A11.getOrder(), u.end(), udummy.begin());
	vdummy = _A22 * udummy;
	copy(vdummy.begin(), vdummy.end(), v.begin() + _A11.getOrder());
	vdummy = _B2 * udummy;
	for (size_t i = 2 * A11.get;i < _B2.numbOfRows;i++)
		v[i] += vdummy[i - (*vdummy.end() - _B2.numbOfRows)];
	// B2^T
	copy(u.end() - _B1.numbOfRows(), u.end(), dummy.begin());
	dummy = B1.t() * dummy;
	for (i = 0; i < A11.getOrder(); ++i)
		v[i] += dummy[i];
}
