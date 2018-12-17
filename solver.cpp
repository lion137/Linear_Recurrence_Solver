// linear recurrence solving using matrix methods

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

#include "solver.h"

using namespace std; 


template <Regular A, Integer N, Operation Op>
A power(A a, N n, Op op) {
	if (n == 0) return identity(a);
		
	while (! odd(n)) {
		a  = op(a, a);
		n = half(n);
	}
	if (n == 1) return a;
	return power_accumulate(a, op(a, a), half(n - 1), op);
}


template <Regular N> 
matrix<N> matr_mul(matrix<N> A, matrix<N> B){ // only square
	int length = A.size();
	matrix<N> C(length, vector<N>(length));
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
			N sum = N(0);
			for (int k = 0; k < length; k++)
				sum = sum +  A[i][k] * B[k][j];
			C[i][j] = sum;
		}
	}
	return C;
}

template <Regular A, Integer N, Operation Op> 
A power_accumulate(A r, A a, N n, Op op) {
	if (n == 0) return r;
	while (true) {
		if (odd(n)) {
			r = op(r, a);
			if (n == 1) return r;
		}
		n  = half(n);
		a = op(a, a);
	}
}

template <Integer N>
bool odd(N n) { return bool (n & 0x1); }

template <MultiplicativeMonoid T> 
matrix<T> identity(matrix<T> e) {
	int length = e.size();
	matrix<T> Out(length, vector<T>(length));
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
			if (i == j) Out[i][j] = 1;
		}
	}
	return Out;
}

template <MultiplicativeMonoid T> 
T identity(T e) {return T(1);}

template <Integer N>
N half (N n) { return n >> 1; }

template <Integer A>
A fibonacci(A n) {
	matrix<A> P = {{A(1), A(1)}, {A(1), A(0)}};
	matrix<A> R = power(P, n, mul_mat<A>());
	return R[1][0];
}

template <Integer N, Integer M> 
N recurrence_solver(M n, std::vector<N> a, std::vector<N> initial){
	int length = a.size();
	matrix<N> A(length, std::vector<N>(length));
	for (int i = 0; i < length; i++)
		A[0][i] = a[i];
	for (int i = 1; i < length; i++) {
		for (int j = 0; j < length; j++) {
			if (i - j == 1)
				A[i][j] = 1;
		}
	}
	reverse(initial.begin(), initial.end());
	return matr_vec_mul(power(A, n + 1 - length, mul_mat<N>()), initial)[0];
}

// matrix * vector, square matrix
template <Regular T>
vector<T> matr_vec_mul(matrix<T> A, vector<T> V) {
	int length = V.size();
	vector<T> Out(length); 
	for (int i = 0; i < length; i++) {
		T sum = T(0);
		for (int k = 0; k < length; k++)
			sum = sum + A[i][k] * V[k];
		Out[i] = sum;
	}
	return Out;
}

// debugging helper, prints a matrix
template <Regular T>void printMatrix(matrix<T> A){
	int length = A.size();
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
			cout << A[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}
	
	
	
	
	
	
	
	
	
