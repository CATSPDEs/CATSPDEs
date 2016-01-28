#include"REAL.h"
#include<vector>
using std::vector;
class Matrix
{
	unsigned int n;							
	vector<REAL> val; //matrix is stored in CRS format
	vector<int> col;
	vector<int> ptr;
private:
	friend float scalar(Matrix&, std::vector<REAL>&);
};
