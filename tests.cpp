// linear recurrence solver tests

#include <iostream>
#include <vector>
#include <cassert>


#include "solver.cpp"
#include "solver.h"

#define run_test(f_name)\
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

long fibonacci2(long n){
	return recurrence_solver(n, vector<long>{1, 1}, vector<long>{0, 1});
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
}

int main () {
	
    cout << "---------------------------\n";
	run_test(testMatrixMul);
	run_test(testPowerNumber);
	run_test(testMatrixPower);
	run_test(testFibonacci);
	run_test(testMatrixVectorMul);
	run_test(testPrintMatrix);
	run_test(testLinearRecurrenceSolver)
	cout << "---------------------------\n";
	cout << "Tests Passed!\n";
	return 0;
}
