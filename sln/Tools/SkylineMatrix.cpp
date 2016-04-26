#include"SkylineMatrix.hpp"

std::istream & SkylineMatrix::loadSparse(std::istream& input)
{
	size_t i, max = _n;
	for (i = 0; i < max; ++i)
		input >> _ptr[i];
	max = _ptr[max];
	for (i = 0; i < max; ++i)
		input >> _col[i];
	for (i = 0; i < max; ++i)
		input >> _diag[i];
	for (i = 0; i < max; ++i)
		input >> _val[i];
	return input;
}

std::ostream& SkylineMatrix::saveSparse(std::ostream& output) const
{
	return output << _n << ' ' << _ptr[_n] + _n << '\n'
		<< _ptr << '\n'
		<< _col << '\n'
		<< _diag << '\n'
		<< _val;
}
