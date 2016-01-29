#include "SLAE.h"
#include<numeric>
#include<cmath>

std::vector<REAL> operator*(const REAL c, const std::vector<REAL>& v)
{
	std::vector<REAL> res(v.size);
	for (int i = 0;i < v.size;i++)
	{
		res[i] = c*v[i];
	}
}

std::vector<REAL> operator+(const std::vector<REAL>& v1, const std::vector<REAL>& v2)
{
	std::vector<REAL> res(v1.size);
	for (int i = 0;i < v1.size;i++)
	{
		res[i] = v1[i]+v2[i];
	}
}

std::vector<REAL> operator-(const std::vector<REAL>& v1, const std::vector<REAL>& v2)
{
	std::vector<REAL> res(v1.size);
	for (int i = 0;i < v1.size;i++)
	{
		res[i] = v1[i] - v2[i];


	}
}

std::vector<REAL> SLAE::CG(REAL e, unsigned iCount)
{
	unsigned n = 0;
	std::vector<REAL> x(this->A.n,0);
	std::vector<REAL> r(this->A.n);
	std::vector<REAL> z;
	REAL a, b, dotProdR, dotProdF;
	r = this->f - A*x;
	z = r;
	dotProdR = std::inner_product(r.begin, r.end, r.begin, 0);
	dotProdF = std::inner_product(this->f.begin, this->f.end, this->f.begin, 0);
	while ((sqrt(dotProdR)/sqrt(dotProdF))>e &&n<iCount)
	{
		a = dotProdR/ std::inner_product((this->A*z).begin, (this->A*z).end, z.begin, .0);
		x = x + a*z;
		r = r - a*(this->A*z);
		b = dotProdR;
		dotProdR= std::inner_product(r.begin, r.end, r.begin, 0);
		b = dotProdR / b;
		z = r + b*z;
		n++;
	}
}

SLAE::SLAE(CRSMatrix & A, std::vector<REAL>& f)
{
	this->A = A;
	this->f = f;
}
