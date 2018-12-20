// linear recurrence solver tests

#include <iostream>
#include <vector>
#include <cassert>
#include <forward_list>

#include "solver.cpp"
#include "solver.h"
//#include "experiments.cpp"
//#include "Matrix.h"
//#include "Matrix.cpp"

#define run_test(f_name);\
	printf("%s\n", #f_name);\
	f_name();

using namespace std;

matrix<int> A = {{1, 1}, {1, 1}};
matrix<int> B = {{1, 1}, {1, 0}};
matrix<int> C = {{2, 1}, {2, 1}};
matrix<int> Id = {{1, 0}, {0, 1}};

void testMatrixMul(){
	assert (C == matr_mul(A, B));
}

void testPowerNumber(){
	assert (power(3, 3, multiplies<int>()) == 27);
	assert (power(1, 0, multiplies<int>()) == 1);
	assert (power(0, 10, multiplies<int>()) == 0);
	assert (power(1, 1, multiplies<int>()) == 1);
	assert (power(34, 1, multiplies<int>()) == 34);
	assert (power(long(10), 11L, multiplies<long>()) == 100000000000L);
}

void testMatrixPower(){
	assert(power(Id, 137, mul_mat<int>()) == identity(Id));	
	
	matrix<int>D = {{4, 4}, {4, 4}};
	assert(power(A, 3, mul_mat<int>()) == D);	
	
	matrix<int> P = {{1, 1}, {1, 0}};
	matrix<int> Z = {{34, 21}, {21, 13}};
	assert(power(P, 8, mul_mat<int>()) == Z);
	
}



void testFibonacci(){
	assert(fibonacci(0) == 0);
	assert(fibonacci(1) ==  1);
	assert(fibonacci(2) ==  1);
	assert(fibonacci(3) ==  2);
	assert(fibonacci(12) ==  144);
	assert(fibonacci(40) ==  102334155);
}


void testMatrixVectorMul(){
	vector<int> V = {0, 0};
	assert (matr_vec_mul(A, V) == V );
	
	vector<int> V1 = {1, 1};
	vector<int> V2 = {2, 2};
	assert( matr_vec_mul(A, V1) == V2);
	matrix<long> A1 = { {1, 2, 3}, {1, 1, 4}, {0, 3, 3}};
	vector<long> V3 = {1, 2, 3};
	vector<long> V4 = {14, 15, 15};
	assert( matr_vec_mul(A1, V3) == V4);
}

void testPrintMatrix(){
	matrix<int> T1(5, vector<int>(5));
	cout << "5 x 5 Zeroes:\n";
	printMatrix(T1);
}

long f(long n){
	if (n == 0) return 0;
	else if (n == 1) return 1;
	else if (n == 2) return 2;
	else 
		return 2 * f(n - 1) + 3 * f(n - 3);
}

 long g(long n){
	if (n == 0) return 0;
	else if (n == 1) return 1;
	else
		return -2 * g(n - 2) - 3 *  g(n - 1) ;
}

void testLinearRecurrenceSolver(){
	vector<long> vv{2, 0, 3}; // recurrence f
	vector<long> in{0, 1, 2};
	assert ( recurrence_solver(30L, vv, in) == f(30L));
	vector<int> vv1{1, 1}; // Fibonacci
	vector<int> in1{0, 1};
	assert (recurrence_solver(40, vv1, in1) == fibonacci(40L));
	vector<long> vv2{-3, -2}; // recurrence g
	vector<long> in2{0, 1};
	assert (recurrence_solver(20, vv2, in2) == g(20));
	vector<long> vv3{3, 0, 0, -5}; // recurrence g
	vector<long> in3{0, 1, 1, 2};
	assert (recurrence_solver(10, vv3, in3) == 1849);
}

void testMulModulo(){
	assert(mul_modulo(4, 10, 7) == 5);
	assert(mul_modulo(4, 10, 5) == 0);
	assert(mul_modulo(0, 0, 5) == 0);
}

template <Operation Op>
int mf(int a, int b, int m, Op op)
	{return op(a, b, m);}

void testMulModuloFunctor(){
	assert(mf(2, 3, 5, mul_mod_functor<int, int>()) == 1);
	assert(mf(7, 7, 7, mul_mod_functor<int, int>()) == 0);
	assert(power_modulo(12345, 123, 100, mul_mod_functor<int, int>()) == 25);
}

void testZeroes(){
	assert (Zero<int>(42) == 0);
	matrix<int> zero = {{0, 0}, {0, 0}};
	assert (Zero<int>(B) == zero);
	
	matrix<unsigned long> C = { {1L, 2L ,3L}, {4L, 5L, 6L}, {0L, 0L, 0L}};
	matrix<unsigned long> ZERO = { {0L, 0L, 0L}, {0L, 0L, 0L}, {0L, 0L, 0L}};
	assert (Zero<unsigned long>(C) == ZERO);
}

void testMulMatrixModulo() {
	matrix<int> A{{2, 4}, {7, 3}};
	matrix<int> B{{1, 2}, {3, 4}};
	matrix<int> C{{4, 0}, {1, 1}};
	assert (mul_mat_modulo(A, B, 5) == C);
	
	matrix<int> AA{{1, 3}, {5, 0}};
	matrix<int> BB{{7, 0}, {2, 3}};
	matrix<int> CC{{1, 1}, {3, 0}};
	assert (mul_mat_modulo(AA, BB, 4) == CC);
}

void testMatrixModulo() {
	matrix<int> A{{1, 2, 30}, {1, 1, 23}, {0, 11, 10}};
	matrix<int> B{{1, 2, 0}, {1, 1, 2}, {0, 2, 1}};
	assert (generic_mod(A, 3) == B);
	assert (generic_mod(42, 3) == 0);
}

void testPowerModulo2() {
	vector<long> a{1, 1};
	vector<long> c{0, 1};
	vector<long> vv{2, 0, 3}; // recurrence f
	vector<long> in{0, 1, 2};
	matrix<int> A{{3, 3}, {2, 1}};
	matrix< long> B{{1L, 1L}, {1L, 0L}};
	for (long i = 1; i < 30; i++)
		assert(power_modulo(B, i - 1, 13L, mul_mat_mod<long, long>())[0][0] == (fibonacci(i) % 13));
	
	matrix<long> C{{2L, 0L, 3L}, {1L, 0L, 0L}, {0L, 1L, 0L}}; // recurrence f matrix (above)
	for (long i = 1; i < 30L; i++){
		assert(power_modulo(C, i - 1, 137L, mul_mat_mod<long, long>())[0][0] == f(i) % 137L);
	}
	
}

long cr(long n) {
	if (n == 0) return 11L;
	else if (n == 1) return 13L;
	else return 5 * cr(n - 1) + 7 * cr(n - 2);
}

void testRecurrenceSolverModulo() {
	vector<long> a{1, 1};
	vector<long> c{0, 1};
	vector<long> vv{2, 0, 3}; // recurrence f
	vector<long> in{0, 1, 2};
	for (long i = 1; i < 40; i++)
		assert(recurrence_solver_modulo(i, a, c, 13L) == (fibonacci(i) % 13));
	
	vector<long> aa{5, 7};
	vector<long> cc{11, 13}; // recurrence cr
	
	for (long i = 0; i < 13; i++){
		assert(recurrence_solver_modulo(i, aa, cc, 13L) == (cr(i) % 13));
	}	
}



int main () { 
    cout << "---------------------------\n";
	run_test(testMatrixMul);
	run_test(testPowerNumber);
	run_test(testMatrixPower);
	run_test(testFibonacci);
	run_test(testMatrixVectorMul);
	run_test(testPrintMatrix);
	run_test(testLinearRecurrenceSolver);
	run_test(testMulModulo);
	run_test(testMulModuloFunctor);
	run_test(testZeroes);
	run_test(testMulMatrixModulo);
	run_test(testMatrixModulo);
	run_test(testPowerModulo2);
	run_test(testRecurrenceSolverModulo);
	cout << "---------------------------\n";
	cout << "Tests Passed!\n";
	return 0;
}
