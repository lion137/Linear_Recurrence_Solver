Solving any linear recurrence relation in ```O(logn)``` time.   
It uses matrix method, as described here: https://en.wikipedia.org/wiki/Recurrence_relation#Solving_via_linear_algebra    
Ex.:    For recurrence:    
```R(n) = 3R(n-1) - 5R(n - 4)```
```R(0) = 0```
```R(1) = 1```
```R(2) = 1```
```R(3) = 2```
Vector ```c``` becomes: ```[3, 0, 0, -5]```, initial vector: ```y = [0, 1, 1, 2]```.    
My ```c``` coefficients vector go from left to right which is is irrelevant. We feed the function ```recurrence solver``` directly.      

Also, solves any linear recurrence modulo m in ```O(logn)``` time. 

Functions are fully generic, so can be extended without problems. I'm using few ```typedef``` shortcuts to easier navigate   
between types, they are not truly Concepts - they are not a part of C++, so far.                    

For longer formulas, like, a hundred coefficients, there is an overhead of matrix multiplication done    
in ```O(n^3)``` time, so it be little bit slower.    

Compiled in gcc c++17, usage(from example above, and Fibonacci):    

```recurrence_solver(10, c, y)``` - returns the 10th value of recurrence above, which is: ```1849```.          

```recurrence_solver_modulo(1234567891011, std::vector<long>{1L, 1L}, std::vector<long>{0L, 1L}, 1000000007L)``` 
gives the ```1234567891011th``` Fibonacci number modulo ```1000000007```, which is ```316399615```. 

Could be used in competitive programming and mathematics, also may help students studying recurrence and/or induction, etc...    







    



