#include"REAL.h"
#include<vector>
using std::vector;
class Matrix
{
	unsigned int n;
	vector<int> diag;
	vector<REAL> upper;
	vector<REAL> lower;
	vector<int> ptr;
	vector<int> rc;

private:
	friend float scalar(Matrix, std::vector<REAL>);
};
