# mfi
Modern Fortran interfaces to BLAS and LAPACK

## TL;DR
```
~$ make test_mfi_blas
~$ ./test_mfi_blas
```

## Support
### BLAS Level 1

|  name  |       description                                        | supported? | 
|--------|----------------------------------------------------------|------------|
| asum   | Sum of vector magnitudes                                 |    [ ]     |
| axpy   | Scalar-vector product                                    |    [ ]     |
| copy   | Copy vector                                              |    [ ]     |
| dot    | Dot product                                              |    [ ]     |
| sdsdot | Dot product with double precision                        |    [ ]     |
| dotc   | Dot product conjugated                                   |    [ ]     |
| dotu   | Dot product unconjugated                                 |    [ ]     |
| nrm2   | Vector 2-norm (Euclidean norm)                           |    [ ]     |
| rot    | Plane rotation of points                                 |    [ ]     |
| rotg   | Generate Givens rotation of points                       |    [ ]     |
| rotm   | Modified Givens plane rotation of points                 |    [ ]     |
| rotmg  | Generate modified Givens plane rotation of points        |    [ ]     |
| scal   | Vector-scalar product                                    |    [ ]     |
| swap   | Vector-vector swap                                       |    [ ]     |
| iamax  | Index of the maximum absolute value element of a vector  |    [ ]     |
| iamin  | Index of the minimum absolute value element of a vector  |    [ ]     |

