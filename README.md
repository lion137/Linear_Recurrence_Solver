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

Functions are fully generic, so can be easily extended, to find recurrence modulo, etc...   
Examples, also, given in ```tests```.        

For longer formulas, like, a hundred coefficients, there is an overhead of matrix multiplication done    
in ```O(n^3)``` time, so it be little bit slower.    



