# DMatrix class
The script matrix.py realizes just the matrix with the operations of multiply, addition, division. There are methods to find determinant, transponent matrix, get zero matrix and etc.

## Properties
+ tp - type of the first element of matrix
+ rows - the amount of rows in matrix
+ columns - the amount of columns
+ shape - the shape of matrix

## Main Methods
+ Init
```
DMatrix(matrix : list) => new DMatrix # initialize matrix 
```
+ Get, set, del item
**To get, set and del items of matrix use simple indexing in squared brackets**
```
m[i][j]
m[i][j][k]
m[j] = val
del m[i][j][k][z]
```
+ Copy
```
copy() => copy of the matrix
```
+ Append
```
append(value) => append value to matrix, and reshape matrix
```
+ Get zero matrix of same shape
```
get_zero(rows = None, columns = None) => returns zero matrix of the same shape if rows = None and columns = None
```
+ Get matrix diagonal
```
diag() => list of diagonal of matrix
```
+ Get the j's row
```
row(j)
```
+ Get the j's column
```
column(j)
```
+ Transpose matrix
```
trans() => returns nothing
```
+ Get the matrix without i's line and j's column
```
mminor(line, row = 0) => DMatrix
```
+ Get the algebraic addition for i's line and j's row
```
alg_adj(i, j=0)
```
+ Get the matrix of algebraic additions of the original matrix
```
malg_adj() => DMatrix
```
+ Get the inverse matrix of original
```
inverse() => DMatrix
```
+ Find determinant of squared matrix
```
det()
```
