#pragma once
#include <vector>

class IRealMatrix { // interface for real-valued matrices that we can do math with
public:
	virtual ~IRealMatrix() {}
	virtual std::vector<double> solve(std::vector<double> const &) = 0; // solve A.x = b
	virtual std::vector<double> mult(std::vector<double> const &) const = 0; // matrix-vector product A.v
};

// useful stuff 
inline std::vector<double> operator*(IRealMatrix const & A, std::vector<double> const & b) { 
	return A.mult(b);
}
inline std::vector<double> operator/(std::vector<double> const & b, IRealMatrix & A) { 
	return A.solve(b);
}
