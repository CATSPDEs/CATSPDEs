// implementation of BandMatrix

#include "BandMatrix.hpp"

BandMatrix::BandMatrix(size_t n, size_t w)
	: AbstractSquareMatrix(n)
	, _w(w)
	, _p((w - 1) / 2)
	, _A(new REAL[n * w])
	, _L(new REAL*[n])
	, _D(&_A[n * _p])
	, _U(new REAL*[n])
	, _isDecomposed(false) {
		if (w < 1) throw std::invalid_argument("width of band must be at least one");
		if (w % 2 == 0) throw std::invalid_argument("width of band cannot be even");
		if (w > 2 * n - 1) throw std::invalid_argument("width of band cannot be greater than number of diagonals");
		unsigned i, l, u;
		for (i = 0, l = 0, u = _n * (_p + 1); i < _n; i++, l += _p, u += _p) {
			_L[i] = &_A[l];
			_U[i] = &_A[u];
		}
}

BandMatrix::~BandMatrix() {
	delete[] _A;
	delete[] _L;
	delete[] _U;
}

// private methods

bool BandMatrix::_computeLU() {
	size_t i, j, k, kMin, m;
	int alj;
	REALs sumL, sumD, sumU;
	/*
		L[i][j] <—> TraditionalMatrix[ali][alj], ali = i,	      alj = i + j - p
		U[i][j] <—> TraditionalMatrix[aui][auj], aui = i + j – p, auj = i 
		alj == aui
	*/
	for (i = 0; i < _n; ++i) {
		sumD = 0.;
		for (j = 0, alj = i - _p, kMin = 0; j < _p; ++j, ++alj) {
			if (alj < 0) {
				++kMin;
				continue;
			}
			sumU = 0.;
			sumL = 0.;
			for (k = kMin, m = _p - j + kMin; k < j; ++k, ++m) {
				sumL += _L[i][k] * _U[alj][m];
				sumU += _U[i][k] * _L[alj][m];
			}
			_L[i][j] -= sumL;
			_L[i][j] /= _D[alj];
			_U[i][j] -= sumU;
			sumD += _L[i][j] * _U[i][j];
		}
		_D[i] -= sumD;
		if (!_D[i]) return false; // we cannot divide by zero
	}
	_isDecomposed = true;
	return true;
}

void BandMatrix::_forwardSubst(std::vector<REAL>& y) const { // L.y = b
	size_t i, j;
	int aj;
	REALs sum;
	// we’re looking through elements of L continuously, from its begining to the last element 
	for (i = 0; i < _n; ++i) {
		sum = 0.;
		for (j = 0, aj = i - _p; j < _p; ++j, ++aj) { // L[i][j] <—> TraditionalMatrix[ai][aj], ai = i, aj = i + j - p
			if (aj < 0) continue;
			sum += _L[i][j] * y[aj];
		}
		y[i] -= sum;
	}
}
/* although approach below is more effective since it has no condition checkings,
we decided to use the one above becuase of clearness and symmetry (look at backSubst method)
REAL* BandMatrix::_forwardSubst(const REAL* b) const {
	unsigned i, shift;
	REAL* y = new REAL[_n];
	y[0] = b[0];
	for (i = 1, shift = _p - 1; shift; ++i, --shift)
		y[i] = b[i] - VectorFunctions::dotProduct(_L[i] + shift, y, i);
	for (shift = 0; i < _n; ++i, ++shift)
		y[i] = b[i] - VectorFunctions::dotProduct(_L[i], y + shift, _p);
	return y;
}
*/

void BandMatrix::_backSubst(std::vector<REAL>& x) const { // U.x = y
	int i, j, ai;
	// we’re looking through elements of U continuously, from its end to the first element 
	for (i = _n - 1; i >= 0; --i) {
		x[i] /= _D[i];
		for (j = _p - 1, ai = i - 1; j >= 0; --j, --ai) { // U[i][j] <—> TraditionalMatrix[ai][aj], ai = i + j – p, aj = i
			if (ai < 0) continue;
			x[ai] -= _U[i][j] * x[i];
		}
	}
}

// public methods

REAL& BandMatrix::set(size_t i, size_t j) {
	if (i >= _n || j >= _n) throw std::out_of_range("matrix doesn’t contain element w/ these indicies");
	if (_diff(i, j) > _p) throw std::invalid_argument("portrait of sparse matrix doesn’t allow to change element w/ these indicies");
	_isDecomposed = false;
	if (i > j) return _L[i][_p - i + j];
	else if (i < j) return _U[j][_p - j + i];
	else return _D[i];
}

REAL BandMatrix::get(size_t i, size_t j) const {
	if (i >= _n || j >= _n) std::out_of_range("matrix doesn’t contain element w/ these indicies");
	if (_diff(i, j) > _p) return 0.;
	if (i > j) return _L[i][_p - i + j];
	else if (i < j) return _U[j][_p - j + i];
	else return _D[i];
}

std::vector<REAL> BandMatrix::solve(std::vector<REAL> const & b) {
	std::vector<REAL> x(_n);
	size_t i;
	if (!_isDecomposed)
		if (!_computeLU()) throw std::runtime_error("LU decomposition doesn’t exist");
	for (i = 0; i < _n; ++i) x[i] = b[i]; // x <— b
	_forwardSubst(x); // x <— L.? = x
	_backSubst(x); // x <— U.? = x
	return x;
}

std::istream& BandMatrix::load(std::istream& input) {
	_isDecomposed = false;
	for (size_t i = 0; i < _n * _w; ++i) input >> _A[i];
	return input;
}

std::ostream& BandMatrix::save(std::ostream& output) const {
	output << _n << ' ' << _w << '\n';
	for (size_t i = 0; i < _n * _w; ++i) output << _A[i] << ' ';
	return output << std::endl;
}

std::vector<REAL> BandMatrix::mult(std::vector<REAL> const & u) const {
	std::vector<REAL> v(_n, 0.);
	// …
	return v;
}
/* TODO: update
REAL* operator*(const BandMatrix& matrix, const REAL* vector) {
	// TODO: check if U & L are null
	unsigned i, j, k;
	REAL* result = new REAL[matrix._n];
	for (i = 0, k = matrix._p; k; ++i, --k) {
		result[i] = matrix._D[i] * vector[i];
		for (j = k; j < matrix._p; ++j) {
			result[i] += matrix._L[i][j] * vector[j - k];
			result[j - k] += matrix._U[i][j] * vector[i];
		}
	}
	for (; i < matrix._n; ++i, ++k) {
		result[i] = matrix._D[i] * vector[i];
		for (j = 0; j < matrix._p; ++j) {
			result[i] += matrix._L[i][j] * vector[j + k];
			result[j + k] += matrix._U[i][j] * vector[i];
		}
	}
	return result;
}
*/