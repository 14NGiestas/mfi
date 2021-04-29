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

| done? | name   | description                                             |
| ----- | ------ | ------------------------------------------------------- |
|       | asum   | Sum of vector magnitudes                                |
| :+1:  | axpy   | Scalar-vector product                                   |
| :+1:  | copy   | Copy vector                                             |
|       | dot    | Dot product                                             |
|       | sdsdot | Dot product with double precision                       |
| :+1:  | dotc   | Dot product conjugated                                  |
|       | dotu   | Dot product unconjugated                                |
|       | nrm2   | Vector 2-norm (Euclidean norm)                          |
|       | rot    | Plane rotation of points                                |
|       | rotg   | Generate Givens rotation of points                      |
|       | rotm   | Modified Givens plane rotation of points                |
|       | rotmg  | Generate modified Givens plane rotation of points       |
|       | scal   | Vector-scalar product                                   |
| :+1:  | swap   | Vector-vector swap                                      |
| :+1:  | iamax  | Index of the maximum absolute value element of a vector |
| :+1:  | iamin  | Index of the minimum absolute value element of a vector |

### BLAS Level 2

| done? | name | description                                                              |
| ----- | ---- | ------------------------------------------------------------------------ |
| :+1:  | gbmv | Matrix-vector product using a general band matrix                        |
| :+1:  | gemv | Matrix-vector product using a general matrix                             |
| :+1:  | ger  | Rank-1 update of a general matrix                                        |
| :+1:  | gerc | Rank-1 update of a conjugated general matrix                             |
| :+1:  | geru | Rank-1 update of a general matrix, unconjugated                          |
| :+1:  | hbmv | Matrix-vector product using a Hermitian band matrix                      |
| :+1:  | hemv | Matrix-vector product using a Hermitian matrix                           |
| :+1:  | her  | Rank-1 update of a Hermitian matrix                                      |
| :+1:  | her2 | Rank-2 update of a Hermitian matrix                                      |
| :+1:  | hpmv | Matrix-vector product using a Hermitian packed matrix                    |
| :+1:  | hpr  | Rank-1 update of a Hermitian packed matrix                               |
| :+1:  | hpr2 | Rank-2 update of a Hermitian packed matrix                               |
| :+1:  | sbmv | Matrix-vector product using symmetric band matrix                        |
| :+1:  | spmv | Matrix-vector product using a symmetric packed matrix                    |
| :+1:  | spr  | Rank-1 update of a symmetric packed matrix                               |
| :+1:  | spr2 | Rank-2 update of a symmetric packed matrix                               |
| :+1:  | symv | Matrix-vector product using a symmetric matrix                           |
| :+1:  | syr  | Rank-1 update of a symmetric matrix                                      |
| :+1:  | syr2 | Rank-2 update of a symmetric matrix                                      |
| :+1:  | tbmv | Matrix-vector product using a triangular band matrix                     |
| :+1:  | tbsv | Solution of a linear system of equations with a triangular band matrix   |
| :+1:  | tpmv | Matrix-vector product using a triangular packed matrix                   |
| :+1:  | tpsv | Solution of a linear system of equations with a triangular packed matrix |
| :+1:  | trmv | Matrix-vector product using a triangular matrix                          |
| :+1:  | trsv | Solution of a linear system of equations with a triangular matrix        |

### BLAS Level 3

| done? | name  | description                                                                                            |
| ----- | ----- | ------------------------------------------------------------------------------------------------------ |
| :+1:  | gemm  | Computes a matrix-matrix product with general matrices.                                                |
| :+1:  | hemm  | Computes a matrix-matrix product where one input matrix is Hermitian and one is general.               |
| :+1:  | herk  | Performs a Hermitian rank-k update.                                                                    |
| :+1:  | her2k | Performs a Hermitian rank-2k update.                                                                   |
| :+1:  | symm  | Computes a matrix-matrix product where one input matrix is symmetric and one matrix is general.        |
| :+1:  | syrk  | Performs a symmetric rank-k update.                                                                    |
| :+1:  | syr2k | Performs a symmetric rank-2k update.                                                                   |
| :+1:  | trmm  | Computes a matrix-matrix product where one input matrix is triangular and one input matrix is general. |
| :+1:  | trsm  | Solves a triangular matrix equation (forward or backward solve).                                       |

### LAPACK

#### Linear Equation Routines

| name  | description                                                                                                                                              |
| ----- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| geqrf | Computes the QR factorization of a general m-by-n matrix.                                                                                                |
| gerqf | Computes the RQ factorization of a general m-by-n matrix.                                                                                                |
| getrf | Computes the LU factorization of a general m-by-n matrix.                                                                                                |
| getri | Computes the inverse of an LU-factored general matrix.                                                                                                   |
| getrs | Solves a system of linear equations with an LU-factored square coefficient matrix, with multiple right-hand sides.                                       |
| hetrf | Computes the Bunch-Kaufman factorization of a complex Hermitian matrix.                                                                                  |
| orgqr | Generates the real orthogonal matrix Q of the QR factorization formed by geqrf.                                                                          |
| ormqr | Multiplies a real matrix by the orthogonal matrix Q of the QR factorization formed by geqrf.                                                             |
| ormrq | Multiplies a real matrix by the orthogonal matrix Q of the RQ factorization formed by gerqf.                                                             |
| potrf | Computes the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix.                                                                 |
| potri | Computes the inverse of a Cholesky-factored symmetric (Hermitian) positive-definite matrix.                                                              |
| potrs | Solves a system of linear equations with a Cholesky-factored symmetric (Hermitian) positive-definite coefficient matrix, with multiple right-hand sides. |
| sytrf | Computes the Bunch-Kaufman factorization of a symmetric matrix.                                                                                          |
| trtrs | Solves a system of linear equations with a triangular coefficient matrix, with multiple right-hand sides.                                                |
| ungqr | Generates the complex unitary matrix Q of the QR factorization formed by geqrf.                                                                          |
| unmqr | Multiplies a complex matrix by the unitary matrix Q of the QR factorization formed by geqrf.                                                             |
| unmrq | Multiplies a complex matrix by the unitary matrix Q of the RQ factorization formed by gerqf.                                                             |

#### Singular Value and Eigenvalue Problem Routines

| name  | description                                                                                                                                             |
| ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| gebrd | Reduces a general matrix to bidiagonal form.                                                                                                            |
| gesvd | Computes the singular value decomposition of a general rectangular matrix.                                                                              |
| heevd | Computes all eigenvalues and, optionally, all eigenvectors of a complex Hermitian matrix using divide and conquer algorithm.                            |
| hegvd | Computes all eigenvalues and, optionally, all eigenvectors of a complex generalized Hermitian definite eigenproblem using divide and conquer algorithm. |
| hetrd | Reduces a complex Hermitian matrix to tridiagonal form.                                                                                                 |
| orgbr | Generates the real orthogonal matrix Q or PT determined by gebrd.                                                                                       |
| orgtr | Generates the real orthogonal matrix Q determined by sytrd.                                                                                             |
| ormtr | Multiplies a real matrix by the orthogonal matrix Q determined by sytrd.                                                                                |
| syevd | Computes all eigenvalues and, optionally, all eigenvectors of a real symmetric matrix using divide and conquer algorithm.                               |
| sygvd | Computes all eigenvalues and, optionally, all eigenvectors of a real generalized symmetric definite eigenproblem using divide and conquer algorithm.    |
| sytrd | Reduces a real symmetric matrix to tridiagonal form.                                                                                                    |
| ungbr | Generates the complex unitary matrix Q or PT determined by gebrd.                                                                                       |
| ungtr | Generates the complex unitary matrix Q determined by hetrd.                                                                                             |
| unmtr | Multiplies a complex matrix by the unitary matrix Q determined by hetrd.                                                                                |
