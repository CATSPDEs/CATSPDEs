#include"BVP.h"

BVP::BVP(REAL k(REAL, REAL), REAL q(REAL, REAL), REAL f(REAL, REAL))
{
	_k = k;
	_q = q;
	_f = f;
}