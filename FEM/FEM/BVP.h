#ifndef BVP_H
#define BVP_H

#include"..\Matrices’ Routine\SLAE.h"

class BVP //BVP stands for boundary value problem
{
	REAL (*_k)(REAL, REAL);
	REAL (*_q)(REAL, REAL);
	REAL (*_f)(REAL, REAL);
	BVP(REAL k(REAL, REAL), REAL q(REAL, REAL), REAL f(REAL, REAL));
};

#endif
