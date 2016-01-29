#include "SLAE.h"
#include <numeric>
#include <cmath>

std::vector<REAL> operator*(const REAL c, const std::vector<REAL>& v)
{
	std::vector<REAL> res(v.size());
	for (int i = 0;i < v.size();i++)
		res[i] = c*v[i];
	return res;
}

std::vector<REAL> operator+(const std::vector<REAL>& v1, const std::vector<REAL>& v2)
{
	std::vector<REAL> res(v1.size());
	for (int i = 0;i < v1.size();i++)
		res[i] = v1[i]+v2[i];
	return res;
}

std::vector<REAL> operator-(const std::vector<REAL>& v1, const std::vector<REAL>& v2)
{
	std::vector<REAL> res(v1.size());
	for (int i = 0;i < v1.size();i++)
		res[i] = v1[i] - v2[i];
	return res;
}

std::vector<REAL> SLAE::_CG(REAL e, unsigned iCount)
{
	unsigned n = 0;
	std::vector<REAL> x(this->_A._n,.0);
	std::vector<REAL> r(this->_A._n);
	std::vector<REAL> z;
	REAL a, b, dotProdR, dotProdF;
	r = this->_f - this->_A*x;
	z = r;
	dotProdR = std::inner_product(r.begin(), r.end(), r.begin(), .0);
	dotProdF = std::inner_product(this->_f.begin(), this->_f.end(), this->_f.begin(), .0);
	while ((sqrt(dotProdR)/sqrt(dotProdF))>e &&n<iCount)
	{
		a = dotProdR/ std::inner_product((this->_A*z).begin(), (this->_A*z).end(), z.begin(), .0);
		x = x + a*z;
		r = r - a*(this->_A*z);
		b = dotProdR;
		dotProdR= std::inner_product(r.begin(), r.end(), r.begin(), .0);
		b = dotProdR / b;
		z = r + b*z;
		n++;
	}
	return x;
}

SLAE::SLAE(CRSMatrix& A, std::vector<REAL>& f) :
	_A(A), _f(f) {
}

std::vector<REAL> SLAE::solve()
{
	unsigned n = 10000000;
	REAL e = 1e-14;
	return _CG(e,n);
}
