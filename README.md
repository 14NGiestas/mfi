# MFI
## Modern Fortran interfaces to BLAS and LAPACK
This project aims to be a collection of modern fortran interfaces to commonly used procedure, for now BLAS and LAPACK.
The main goal is to reduce the pain of using such libraries, providing a generic interface to the intrinsic supported types, 
identifying the optional or reconstructible arguments of a given procedure. The code uses [fypp](https://github.com/aradi/fypp),
to generate the interfaces automatically to all supported types and kinds.

## TL;DR
```
$ git clone https://github.com/14NGiestas/mfi.git
$ sudo pip install fypp
$ cd mfi/
$ make test_mfi_blas
$ ./test_mfi_blas
```

## Support
### BLAS Level 1

|  name  |       description                                        | done? | 
|--------|----------------------------------------------------------|-------|
| asum   | Sum of vector magnitudes                                 |       |
| axpy   | Scalar-vector product                                    |  :+1: |
| copy   | Copy vector                                              |       |
| dot    | Dot product                                              |       |
| sdsdot | Dot product with double precision                        |       |
| dotc   | Dot product conjugated                                   |  :+1: |
| dotu   | Dot product unconjugated                                 |       |
| nrm2   | Vector 2-norm (Euclidean norm)                           |       |
| rot    | Plane rotation of points                                 |       |
| rotg   | Generate Givens rotation of points                       |       |
| rotm   | Modified Givens plane rotation of points                 |       |
| rotmg  | Generate modified Givens plane rotation of points        |       |
| scal   | Vector-scalar product                                    |       |
| swap   | Vector-vector swap                                       |       |
| iamax  | Index of the maximum absolute value element of a vector  |  :+1: |
| iamin  | Index of the minimum absolute value element of a vector  |  :+1: |

### BLAS Level 2

|  name  |       description                                                         | done? | 
|--------|---------------------------------------------------------------------------|-------|
|  gbmv  | Matrix-vector product using a general band matrix                         |       |
|  gemv  | Matrix-vector product using a general matrix                              |  :+1: |
|  ger   | Rank-1 update of a general matrix                                         |       |
|  gerc  | Rank-1 update of a conjugated general matrix                              |       |
|  geru  | Rank-1 update of a general matrix, unconjugated                           |       |
|  hbmv  | Matrix-vector product using a Hermitian band matrix                       |       |
|  hemv  | Matrix-vector product using a Hermitian matrix                            |       |
|  her   | Rank-1 update of a Hermitian matrix                                       |       |
|  her2  | Rank-2 update of a Hermitian matrix                                       |       |
|  hpmv  | Matrix-vector product using a Hermitian packed matrix                     |       |
|  hpr   | Rank-1 update of a Hermitian packed matrix                                |       |
|  hpr2  | Rank-2 update of a Hermitian packed matrix                                |       |
|  sbmv  | Matrix-vector product using symmetric band matrix                         |       |
|  spmv  | Matrix-vector product using a symmetric packed matrix                     |       |
|  spr   | Rank-1 update of a symmetric packed matrix                                |       |
|  spr2  | Rank-2 update of a symmetric packed matrix                                |       |
|  symv  | Matrix-vector product using a symmetric matrix                            |       |
|  syr   | Rank-1 update of a symmetric matrix                                       |       |
|  syr2  | Rank-2 update of a symmetric matrix                                       |       |
|  tbmv  | Matrix-vector product using a triangular band matrix                      |       |
|  tbsv  | Solution of a linear system of equations with a triangular band matrix    |       |
|  tpmv  | Matrix-vector product using a triangular packed matrix                    |       |
|  tpsv  | Solution of a linear system of equations with a triangular packed matrix  |       |
|  trmv  | Matrix-vector product using a triangular matrix                           |       |
|  trsv  | Solution of a linear system of equations with a triangular matrix         |       |


### BLAS Level 3

|  name  |       description                                                                                      | done? | 
|--------|--------------------------------------------------------------------------------------------------------|-------|
|  gemm  | Computes a matrix-matrix product with general matrices.                                                |  :+1: |
|  hemm  | Computes a matrix-matrix product where one input matrix is Hermitian and one is general.               |       |
|  herk  | Performs a Hermitian rank-k update.                                                                    |  :+1: |
| her2k  | Performs a Hermitian rank-2k update.                                                                   |       |
|  symm  | Computes a matrix-matrix product where one input matrix is symmetric and one matrix is general.        |       |
|  syrk  | Performs a symmetric rank-k update.                                                                    |       |
| syr2k  | Performs a symmetric rank-2k update.                                                                   |       |
|  trmm  | Computes a matrix-matrix product where one input matrix is triangular and one input matrix is general. |       |
|  trsm  | Solves a triangular matrix equation (forward or backward solve).                                       |       |
