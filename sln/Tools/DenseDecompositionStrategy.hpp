#pragma once
#include "vector.hpp"

namespace DenseDecompositionStrategy {

	// row-oriented IKJ (L + I) . (U + D) decomposition from Saad, p. 305
	template <typename T>
	inline void IKJ(std::vector<std::vector<T>>& A) {
		Index n = A.size();
		for (Index i = 1; i < n; ++i)
			for (Index k = 0; k < i; ++k) {
				A[i][k] /= A[k][k];
				for (Index j = k + 1; j < n; ++j)
					A[i][j] -= A[i][k] * A[k][j];
			}
	}

	// Crout (L + I) . (U + D) decomposition from Saad, p. 332
	template <typename T>
	inline void Crout(std::vector<std::vector<T>>& A) {
		Index n = A.size();
		for (Index k = 0; k < n; ++k) {
			for (Index i = 0; i < k; ++i) {
				if (A[k][i]) for (Index j = k; j < n; ++j) A[k][j] -= A[k][i] * A[i][j];
				if (A[i][k]) for (Index j = k + 1; j < n; ++j) A[j][k] -= A[i][k] * A[j][i];
			}
			for (Index i = k + 1; i < n; ++i) A[i][k] /= A[k][k];
		}
	}

	template <typename T>
	inline void Crout_L_plus_D_times_D_plus_U(std::vector<std::vector<T>>& A) {
		Index n = A.size();
		for (Index k = 0; k < n; ++k) {
			for (Index i = 0; i < k; ++i) {
				A[k][k] -= A[k][i] * A[i][k];
				if (A[k][i]) for (Index j = k + 1; j < n; ++j) A[k][j] -= A[k][i] * A[i][j];
				if (A[i][k]) for (Index j = k + 1; j < n; ++j) A[j][k] -= A[i][k] * A[j][i];
			}
			A[k][k] = sqrt(A[k][k]);
			for (Index i = k + 1; i < n; ++i) {
				A[i][k] /= A[k][k];
				A[k][i] /= A[k][k];
			}
		}
	}
	
	template <typename T>
	inline void Crout_L_plus_I_times_D_times_I_plus_U(std::vector<std::vector<T>>& A) {
		Index n = A.size();
		for (Index k = 0; k < n; ++k) {
			for (Index i = 0; i < k; ++i) {
				A[k][k] -= A[i][i] * A[k][i] * A[i][k];
				if (A[k][i]) for (Index j = k + 1; j < n; ++j) A[k][j] -= A[i][i] * A[k][i] * A[i][j];
				if (A[i][k]) for (Index j = k + 1; j < n; ++j) A[j][k] -= A[i][i] * A[i][k] * A[j][i];
			}
			for (Index i = k + 1; i < n; ++i) {
				A[i][k] /= A[k][k];
				A[k][i] /= A[k][k];
			}
		}
	}

}
