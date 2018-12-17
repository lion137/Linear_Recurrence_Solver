// header

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

#define Integer typename
#define Operation typename
#define Regular typename
#define MultiplicativeMonoid typename	



template <Regular A> using matrix = std::vector<std::vector<A> >;

template <Regular N> matrix<N> matr_mul(matrix<N> A, matrix<N> B); // matrix multiplication

template <Regular A, Integer N, Operation Op> 
A power_accumulate(A r, A a, N n, Op op); // accumulator - power helper 

template <Regular T> struct mul_mat { //matrix multilication functor
  matrix<T> operator() (const matrix<T>& x, const matrix<T>& y)  {return matr_mul(x, y);}
 };

template <Integer N>
bool odd(N n);

template <Integer N>
N half (N n); 

template <MultiplicativeMonoid T> 
T identity(T e);

template <MultiplicativeMonoid T> 
matrix<T> identity(matrix<T> e);

template <Regular A, Integer N, Operation Op>
A power(A a, N n, Op op);

template <Integer A> A fibonacci(A n);

// linear recurrence solver, takes an integer n, a vector of parameters,
// a vector of initial values and returns the nth value of a relation
template <Integer N, Integer M> N recurrence_solver(M n, std::vector<N> a, std::vector<N> initial);


// matrix * vector
template <Regular T>
std::vector<T> matr_vec_mul(matrix<T> A, std::vector<T> V);

// debugging helper, prints matrix
template <Regular T>void printMatrix(matrix<T> A);

#endif
