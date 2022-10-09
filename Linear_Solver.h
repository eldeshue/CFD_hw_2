#ifndef _LINEAR_SOLVER_
#define _LINEAR_SOLVER_

#include "Dense_Matrix.h"

// Gauss Elimination

// LU Decomposition

// Cholesky Algorithm

// Tohmas Algorithm
// Solving A*X = B, with Tridiagonal Matrix A
template<typename T, int size>
Dense_Matrix<T, size, 1> Thomas_Algorithm( Dense_Matrix<T, size, size> A, Dense_Matrix<T, size, 1> B) {// call by value
	T result[size][1];
	// step 1. Modifiying Matrix
	// divide by diagonal term (i, i)
	// delete (i + 1, i) term by previous calculation
	// repeat this from 0 to size - 1.
	for (int i = 0; i < size - 1; i++) {
		// divide by diagonal term A.data_table[i][i]
		B.data_table[i][0] /= A.data_table[i][i];
		A.data_table[i][i + 1] /= A.data_table[i][i];
		A.data_table[i][i] = 1;
		// delete (i+1, i term) term
		B.data_table[i + 1][0] -= A.data_table[i + 1][i] * B.data_table[i][0];
		A.data_table[i + 1][i + 1] -= A.data_table[i + 1][i] * A.data_table[i][i + 1];
		A.data_table[i + 1][i] = 0;
	}
	B.data_table[size - 1][0] /= A.data_table[size - 1][size - 1];
	A.data_table[size - 1][size - 1] = 1;
	//A.show();

	// step 2. Result Calculation
	// from size - 1 to 0, calculate the value of result matrix (i, 0).
	// 이후 이를 다음 행에 반영시켜줌.
	result[size - 1][0] = B.data_table[size - 1][0];
	for (int i = size - 2; i >= 0; i--) {
		result[i][0] = ( B.data_table[i][0] 
			- A.data_table[i][i + 1] * result[i + 1][0] );
	}
	// return result
	return Dense_Matrix<T, size, 1>(result);
}

#endif
