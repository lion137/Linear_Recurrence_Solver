// linear recurrence solving using matrix methods

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

#include "solver.h"

using namespace std; 


template <Regular A, Integer N, Operation Op> // power
A power(A a, N n, Op op) {
	if (n == 0) return identity(a);
		
	while (! odd(n)) {
		a  = op(a, a);
		n = half(n);
	}
	if (n == 1) return a;
	return power_accumulate(a, op(a, a), half(n - 1), op);
}

template <Integer N, Integer M>  //recurrence solver
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

template <Integer N>  //recurrence solver modulo
N recurrence_solver_modulo(N n, std::vector<N> a, std::vector<N> initial, N mod){
	N length(a.size());
	matrix<N> A(length, std::vector<N>(length));
	for (int i = 0; i < length; i++)
		A[0][i] = a[i];
	for (int i = 1; i < length; i++) {
		for (int j = 0; j < length; j++) {
			if (i - j == 1)
				A[i][j] = 1;
		}
	}
	
	if ( n < initial.size()) return initial[n] % mod;
	reverse(initial.begin(), initial.end());
	return matr_vec_mul_mod(power_modulo(A, n + 1 - length, mod,  mul_mat_mod<N, N>()), initial, mod)[0];
}


template <Regular A, Integer N, Operation Op> // power modulo 
A power_modulo(A a, N n, N mod, Op op) {
	if (mod == 1) return Zero(a);
	if (a == Zero(a)) return Zero(a);
	A result(identity(a));
	a = generic_mod(a, mod);
	while (n > 0) {
		if (odd(n)) {
			result = op(a, result, mod);
		}
		n = half(n);
		a = op(a, a, mod);
	}
	return result;
}

template <Integer A>
A fibonacci(A n) {
	matrix<A> P = {{A(1), A(1)}, {A(1), A(0)}};
	matrix<A> R = power(P, n, mul_mat<A>());
	return R[1][0];
}

/* --------------- HELPER FUNCTIONS ----------------------*/

template <Regular T>
T Zero(T e) {
	return T(0);
}

template <RegularMatrix T>
matrix<T> Zero(matrix<T> e) {
	return matrix<T>(e.size(), vector<T>(e.size()));
}

template <Regular T, Integer N>// possible overflow
T mul_modulo(T a, T b, N mod){
	return ( (a % mod) * (b % mod)) % mod;
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

//matrix modulo multiplication
template <Regular T, Integer N> matrix<T> mul_mat_modulo(matrix<T> A, matrix<T> B, N mod){
	int length = A.size();
	matrix<T> C(length, vector<T>(length));
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
			N sum = N(0);
			for (int k = 0; k < length; k++)
				sum = (sum +  mul_modulo(A[i][k], B[k][j], mod)) % mod;// possible overflow
			C[i][j] = sum;
		}
	}
	return C;
}

template <Regular T, Integer N> // matrix mod
matrix<T> generic_mod(matrix<T> A, N m){
	N length = A.size();
	for (int i = 0; i < length; i++){
		for (int j = 0; j < length; j++){
				A[i][j] = A[i][j] % m;
		}
	}
	return A;
}

template <Regular T, Integer N>
T generic_mod(T t, N m) {
	return t % m;
}

template <Regular A, Integer N, Operation Op> 
A power_accumulate(A r, A a, N n, Op op) {
	if (n == 0) return r;
	while (true) {
		if (odd(n)) {
			r = op(a, r);
			if (n == 1) return r;
		}
		n  = half(n);
		a = op(a, a);
	}
}

template <Regular A, Integer N, Operation Op> 
A power_accumulate_mod(A r, A a, N n, N mod, Op op) {
	if (n == 0) return r;
	while (true) {
		if (odd(n)) {
			r = op(r, a, mod);
			if (n == 1) return r;
		}
		n  = half(n);
		a = op(a, a, mod);
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

template <Regular T, Integer N> // matr * vec mod 
vector<T> matr_vec_mul_mod(matrix<T> A, vector<T> V, N m) {
	int length = V.size();
	vector<T> Out(length); 
	for (int i = 0; i < length; i++) {
		T sum = T(0);
		for (int k = 0; k < length; k++)
			sum = generic_mod(generic_mod(sum, m) + mul_modulo(A[i][k], V[k], m), m);
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
	
	

	
	
	
	
	
