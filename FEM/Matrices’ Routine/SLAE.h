#ifndef SLAE_H
#define SLAE_H

#include"CRSMatrix.h"

class SLAE
{
	CRSMatrix A;
	std::vector<REAL> f;
	std::vector<REAL> CG(REAL, unsigned); //conjugate gradients method
public:
	SLAE(CRSMatrix&, std::vector<REAL>&);

};
#endif