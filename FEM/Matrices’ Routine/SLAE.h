#ifndef SLAE_H
#define SLAE_H

#include "CRSMatrix.h"

class SLAE
{
	CRSMatrix _A;
	std::vector<REAL> _f;
	std::vector<REAL> _CG(REAL, unsigned); //conjugate gradients method
public:
	SLAE(CRSMatrix&, std::vector<REAL>&);
	std::vector<REAL> solve();
};
#endif