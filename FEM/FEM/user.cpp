#include"BVP.h"
#include<iostream>
REAL k_(REAL x, REAL y)
{
	return 1;
}

REAL q_(REAL x, REAL y)
{
	return 2;
}

REAL f_(REAL x, REAL y)
{
	return 3;
}

int main()
{
	REAL(*k)(REAL x, REAL y);
	REAL(*q)(REAL x, REAL y);
	REAL(*f)(REAL x, REAL y);
	BVP problem(k, q, f);
}