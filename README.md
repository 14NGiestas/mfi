# MFI

## Modern Fortran interfaces to BLAS and LAPACK

This project aims to be a collection of modern fortran interfaces to commonly used procedure, for now BLAS and LAPACK.
The main goal is to reduce the pain of using such libraries, providing a generic interface to the intrinsic supported types and 
identifying the optional or reconstructible arguments of a given procedure. The code uses [fypp](https://github.com/aradi/fypp),
to generate the interfaces automatically to all supported types and kinds.

### Example $C = AB$

```fortran
program main
use mfi_blas, only: mfi_gemm
use f77_blas, only: f77_gemm
use iso_fortran_env
implicit none
! ... variables and other boilerplate code here ...
! Original interface: type and precision dependent
call cgemm('N','N', N, N, N, alpha, A, N, B, N, beta, C, N)
! Improved F77 interface: still a lot of arguments
call f77_gemm('N','N', N, N, N, alpha, A, N, B, N, beta, C, N)
! Modern fortran interface: less arguments and more readable 
call mfi_gemm(A,B,C)
end program
```

## Getting Started

First get the code, by cloning the repo:

```sh
git clone https://github.com/14NGiestas/mfi.git
cd mfi/
```

Install the [fypp](https://github.com/aradi/fypp) using the command:

```sh
sudo pip install fypp
```

### FPM

This project supports the [Fortran Package Manager](https://github.com/fortran-lang/fpm). Follow the directions on that page to install FPM if you haven't already.

### Running 
```sh
make
fpm test
```

## Support

### BLAS Level 1

| name   | description                                             | done? |
| ------ | ------------------------------------------------------- | ----- |
| asum   | Sum of vector magnitudes                                |       |
| axpy   | Scalar-vector product                                   | :+1:  |
| copy   | Copy vector                                             | :+1:  |
| dot    | Dot product                                             |       |
| sdsdot | Dot product with double precision                       |       |
| dotc   | Dot product conjugated                                  | :+1:  |
| dotu   | Dot product unconjugated                                |       |
| nrm2   | Vector 2-norm (Euclidean norm)                          |       |
| rot    | Plane rotation of points                                |       |
| rotg   | Generate Givens rotation of points                      |       |
| rotm   | Modified Givens plane rotation of points                |       |
| rotmg  | Generate modified Givens plane rotation of points       |       |
| scal   | Vector-scalar product                                   |       |
| swap   | Vector-vector swap                                      | :+1:  |
| iamax  | Index of the maximum absolute value element of a vector | :+1:  |
| iamin  | Index of the minimum absolute value element of a vector | :+1:  |

### BLAS Level 2

| name | description                                                              | done? |
| ---- | ------------------------------------------------------------------------ | ----- |
| gbmv | Matrix-vector product using a general band matrix                        | :+1:  |
| gemv | Matrix-vector product using a general matrix                             | :+1:  |
| ger  | Rank-1 update of a general matrix                                        | :+1:  |
| gerc | Rank-1 update of a conjugated general matrix                             | :+1:  |
| geru | Rank-1 update of a general matrix, unconjugated                          | :+1:  |
| hbmv | Matrix-vector product using a Hermitian band matrix                      | :+1:  |
| hemv | Matrix-vector product using a Hermitian matrix                           | :+1:  |
| her  | Rank-1 update of a Hermitian matrix                                      | :+1:  |
| her2 | Rank-2 update of a Hermitian matrix                                      | :+1:  |
| hpmv | Matrix-vector product using a Hermitian packed matrix                    |       |
| hpr  | Rank-1 update of a Hermitian packed matrix                               |       |
| hpr2 | Rank-2 update of a Hermitian packed matrix                               |       |
| sbmv | Matrix-vector product using symmetric band matrix                        | :+1:  |
| spmv | Matrix-vector product using a symmetric packed matrix                    |       |
| spr  | Rank-1 update of a symmetric packed matrix                               |       |
| spr2 | Rank-2 update of a symmetric packed matrix                               |       |
| symv | Matrix-vector product using a symmetric matrix                           | :+1:  |
| syr  | Rank-1 update of a symmetric matrix                                      | :+1:  |
| syr2 | Rank-2 update of a symmetric matrix                                      | :+1:  |
| tbmv | Matrix-vector product using a triangular band matrix                     | :+1:  |
| tbsv | Solution of a linear system of equations with a triangular band matrix   | :+1:  |
| tpmv | Matrix-vector product using a triangular packed matrix                   |       |
| tpsv | Solution of a linear system of equations with a triangular packed matrix |       |
| trmv | Matrix-vector product using a triangular matrix                          | :+1:  |
| trsv | Solution of a linear system of equations with a triangular matrix        | :+1:  |

### BLAS Level 3

| name  | description                                                                                            | done? |
| ----- | ------------------------------------------------------------------------------------------------------ | ----- |
| gemm  | Computes a matrix-matrix product with general matrices.                                                | :+1:  |
| hemm  | Computes a matrix-matrix product where one input matrix is Hermitian and one is general.               | :+1:  |
| herk  | Performs a Hermitian rank-k update.                                                                    | :+1:  |
| her2k | Performs a Hermitian rank-2k update.                                                                   | :+1:  |
| symm  | Computes a matrix-matrix product where one input matrix is symmetric and one matrix is general.        | :+1:  |
| syrk  | Performs a symmetric rank-k update.                                                                    | :+1:  |
| syr2k | Performs a symmetric rank-2k update.                                                                   | :+1:  |
| trmm  | Computes a matrix-matrix product where one input matrix is triangular and one input matrix is general. | :+1:  |
| trsm  | Solves a triangular matrix equation (forward or backward solve).                                       | :+1:  |
