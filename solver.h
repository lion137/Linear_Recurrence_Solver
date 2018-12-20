// header

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

#define Integer typename
#define Operation typename
#define Regular typename
#define RegularMatrix typename
#define MultiplicativeMonoid typename	



template <Regular A> using matrix = std::vector<std::vector<A> >;

template <Regular A, Integer N, Operation Op> // power
A power(A a, N n, Op op);

// linear recurrence solver, takes an integer n, a vector of parameters,
// a vector of initial values and returns the nth value of a relation
template <Integer N, Integer M> N recurrence_solver(M n, std::vector<N> a, std::vector<N> initial);

// linear recurrence solver modulo
template <Integer N>  //recurrence solver
N recurrence_solver_modulo(N n, std::vector<N> a, std::vector<N> initial, N mod);

template <Regular A, Integer N, Operation Op> // power modulo 
A power_modulo(A a, N n, N mod, Op op);

/* ------------ MATRIX FUNCTIONS ----------------*/

template <Regular N> matrix<N> matr_mul(matrix<N> A, matrix<N> B); // matrix multiplication

template <Regular T> struct mul_mat { //matrix multilication functor
  matrix<T> operator() (const matrix<T>& x, const matrix<T>& y)  {return matr_mul(x, y);}
 };

template <Regular T, Integer N> matrix<T> mul_mat_modulo(matrix<T> A, matrix<T> B, N mod);

template <Regular T, Integer N> struct mul_mat_mod { //matrix multilication functor modulo
  matrix<T> operator() (const matrix<T>& x, const matrix<T>& y, const N& m)  {return mul_mat_modulo(x, y, m);}
 };

template <Regular T, Integer N>
T mul_modulo(T a, T b, N mod);

template <Regular T, Integer N> struct mul_mod_functor { //mult modulo functor
  T operator() (const T& x, const T& y, const N& m)  { return mul_modulo(x, y, m);}
};

// matrix * vector
template <Regular T>
std::vector<T> matr_vec_mul(matrix<T> A, std::vector<T> V);

// matrix * vec mod
template <Regular T, Integer N>
std::vector<T> matr_vec_mul_mod(matrix<T> A, std::vector<T> V, N m);

template <Regular T, Integer N>
matrix<T> generic_mod(matrix<T> A, N m);

/*-------------------END MATRIX -------------------*/



/* ------------------HELPERS ----------------------*/

template <Regular T, Integer N>
T generic_mod(T t, N m);

template <Regular T, Integer N>
T mul_modulo(T a, T b, N mod);

template <Regular T>
T Zero(T e);

template <RegularMatrix T>
matrix<T> Zero(matrix<T> e);

template <Integer N>
bool odd(N n);

template <Integer N>
N half (N n); 

template <MultiplicativeMonoid T> 
T identity(T e);

template <MultiplicativeMonoid T> 
matrix<T> identity(matrix<T> e);

template <Integer A> A fibonacci(A n);

template <Regular A, Integer N, Operation Op> 
A power_accumulate(A r, A a, N n, Op op); // accumulator - power helper

template <Regular A, Integer N, Operation Op> // modulo accumulator
A power_accumulate_mod(A r, A a, N n, N mod, Op op);

// debugging helper, prints matrix
template <Regular T>void printMatrix(matrix<T> A);

#endif
