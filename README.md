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

### FPM

This project supports the [Fortran Package Manager](https://github.com/fortran-lang/fpm).
Follow the directions on that page to install FPM if you haven't already.

### Using as a dependency in FPM

Add a entry in the "dependencies" section of your project's fpm.toml

```toml
# fpm.toml
[ dependencies ]
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-fpm" }
```

### Manual building

First get the code, by cloning the repo:

```sh
git clone https://github.com/14NGiestas/mfi.git
cd mfi/
```

### Dependencies

Install the [fypp](https://github.com/aradi/fypp) using the command:

```sh
sudo pip install fypp
```

Install lapack and blas (use the static versions).
This can be tricky, if you run into any problem, please open an issue.

#### Arch Linux
- [aur/openblas-lapack-static](https://aur.archlinux.org/packages/openblas-lapack-static)

#### Ubuntu
- [lapack-dev](https://packages.ubuntu.com/search?suite=default&section=all&arch=any&keywords=lapack-dev&searchon=names)

Usually you can do the following:

```sh
make
fpm test
```

By default, the `lapack-dev` package (which provides the reference blas) do not provide the `i?amin` implementation (among other extensions)
in such cases you can use blas extensions with:

```sh
make FYPPFLAGS=-DMFI_EXTENSIONS
fpm test
```

Or if you have support to such extensions in your blas provider you can:

```sh
make FYPPFLAGS="-DMFI_EXTENSIONS -DMFI_LINK_EXTERNAL"
fpm test
```

which will generate the code linking extensions to the external library


## Support

Please note that this project is experimental, errors and mistakes are to be expected.

There four levels of interfaces that can be used:

1. original f77: explicit declared original interface.
```fortran
call cgemm('N','N', N, N, N, alpha, A, N, B, N, beta, C, N)
```
2. improved f77: original argument convention without need of a prefix.
```fortran
call f77_gemm('N','N', N, N, N, alpha, A, N, B, N, beta, C, N)
```
3. modern interface with prefix:
```fortran
call mfi_sgemm(A,B,C)
```
4. modern interface:
```fortran
call mfi_gemm(A,B,C)
```

If you are searching for a specific interface check the [API reference](https://14ngiestas.github.io/mfi/)



### BLAS
#### Level 1
Most of BLAS level 1 routines can be replaced by intrinsincs and other features in modern fortran.
<details>

|done| name   | description                                             | modern alternative |
|----| ------ | ------------------------------------------------------- | ------------------ |
|:+1:| asum   | Sum of vector magnitudes                                | [sum](https://gcc.gnu.org/onlinedocs/gfortran/SUM.html) |
|:+1:| axpy   | Scalar-vector product                                   | `a*x + b` |
|:+1:| copy   | Copy vector                                             |  `x = b`  |
|:+1:| dot    | Dot product                                             | [dot_product](https://gcc.gnu.org/onlinedocs/gfortran/DOT_005fPRODUCT.html)   |
|:+1:| dotc   | Dot product conjugated                                  | |
|:+1:| dotu   | Dot product unconjugated                                | |
|og77| sdsdot | Compute the inner product of two vectors with extended precision accumulation.            | |
|og77| dsdot  | Compute the inner product of two vectors with extended precision accumulation and result. | |
|:+1:| nrm2   | Vector 2-norm (Euclidean norm)                          | [norm2](https://gcc.gnu.org/onlinedocs/gfortran/NORM2.html) |
|:+1:| rot    | Plane rotation of points                                | |
|:+1:| rotg   | Generate Givens rotation of points                      | |
|:+1:| rotm   | Modified Givens plane rotation of points                | |
|:+1:| rotmg  | Generate modified Givens plane rotation of points       | |
|:+1:| scal   | Vector-scalar product                                   | `a*x + b` |
|:+1:| swap   | Vector-vector swap                                      | |
</details>

#### Level 1 - Utils / Extensions
<details>

| done? | name  | description                                              |  modern alternatives | obs |
| ----- | ----- | -------------------------------------------------------- | ------------------- | --- |
| :+1:  | iamax | Index of the maximum absolute value element of a vector  | [maxval](https://gcc.gnu.org/onlinedocs/gfortran/MAXVAL.html), [maxloc](https://gcc.gnu.org/onlinedocs/gfortran/MAXLOC.html) | |
| :+1:  | iamin | Index of the minimum absolute value element of a vector  | [minval](https://gcc.gnu.org/onlinedocs/gfortran/MINVAL.html), [minloc](https://gcc.gnu.org/onlinedocs/gfortran/MINLOC.html) | |
| :+1:  | lamch | Determines precision machine parameters.                 | [huge](https://gcc.gnu.org/onlinedocs/gfortran/intrinsic-procedures/huge.html), [tiny](https://gcc.gnu.org/onlinedocs/gfortran/intrinsic-procedures/tiny.html), [epsilon](https://gcc.gnu.org/onlinedocs/gfortran/intrinsic-procedures/epsilon.html) | Obs: had to add a parameter so fortran can distinguish between the single and double precision with the same interface. For values of cmach see: [lamch](https://www.netlib.org/lapack//explore-html/d4/d86/group__lamch.html)|
</details>

#### Level 2

<details>

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
</details>

#### Level 3

<details>

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

</details>

### LAPACK :warning:

- Lapack is really huge, so I'm going to focus on getting the improved f77 interfaces ready first.
  Anything I end up using I'm going to implement.

#### Linear solve, $AX = B$
<details>
<!-- ##### LU: General matrix, driver -->

<!-- ##### LU: computational routines (factor, cond, etc.) -->
 
<!-- ##### Cholesky: Hermitian/symmetric positive definite matrix, driver -->

##### Cholesky: computational routines (factor, cond, etc.)
| done?| name  | description               |
| ---- | ----- | ------------------------- |
| :+1: | pocon | condition number estimate |

<!-- ##### LDL: Hermitian/symmetric indefinite matrix, driver -->
 
<!-- ##### LDL: computational routines (factor, cond, etc.) -->
 
<!-- ##### Triangular computational routines (solve, cond, etc.) -->
 
<!-- ##### Auxiliary routines -->
</details>

##### Orthogonal/unitary factors (QR, CS, etc.)
<details>

| done? | name  | description  |
| ----- | ----- | ------------ |
| :+1:  | geqrf | Computes the QR factorization of a general m-by-n matrix. |
| :+1:  | gerqf | Computes the RQ factorization of a general m-by-n matrix. |
| :+1:  | getrf | Computes the LU factorization of a general m-by-n matrix. |
| :+1:  | getri | Computes the inverse of an LU-factored general matrix.    |
| :+1:  | getrs | Solves a system of linear equations with an LU-factored square coefficient matrix, with multiple right-hand sides. |
| :+1:  | hetrf | Computes the Bunch-Kaufman factorization of a complex Hermitian matrix. |
| :+1:  | potrf | Computes the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix.     |
| :+1:  | potri | Computes the inverse of a Cholesky-factored symmetric (Hermitian) positive-definite matrix.  |
| :+1:  | potrs | Solves a system of linear equations with a Cholesky-factored symmetric (Hermitian) positive-definite coefficient matrix, with multiple right-hand sides.  |
| :+1:  | orgqr | Generates the real orthogonal matrix Q of the QR factorization formed by geqrf. |
| :+1:  | orgrq | Generates the real orthogonal matrix Q of the RQ factorization formed by gerqf. |
| :+1:  | ormqr | Multiplies a real matrix by the orthogonal matrix Q of the QR factorization formed by geqrf. |
| f77  | ormrq | Multiplies a real matrix by the orthogonal matrix Q of the RQ factorization formed by gerqf. |
| :+1:  | sytrf | Computes the Bunch-Kaufman factorization of a symmetric matrix.                        |
| :+1:  | trtrs | Solves a system of linear equations with a triangular coefficient matrix, with multiple right-hand sides. |
| :+1:  | ungqr | Generates the complex unitary matrix Q of the QR factorization formed by geqrf.  |
| :+1:  | ungrq | Generates the complex unitary matrix Q of the RQ factorization formed by gerqf.  |
| :+1:  | unmqr | Multiplies a complex matrix by the unitary matrix Q of the QR factorization formed by geqrf. |
| f77  | unmrq | Multiplies a complex matrix by the unitary matrix Q of the RQ factorization formed by gerqf. |
| :+1:  | org2r | Generates the real orthogonal matrix Q of the QR factorization formed by geqr2. |
| :+1:  | orm2r | Multiplies a real matrix by the orthogonal matrix Q formed by geqr2. |
| :+1:  | ung2r | Generates the complex unitary matrix Q of the QR factorization formed by geqr2. |
| :+1:  | unm2r | Multiplies a complex matrix by the unitary matrix Q formed by geqr2. |
| :+1:  | orgr2 | Generates the real orthogonal matrix Q of the RQ factorization formed by gerq2. |
| :+1:  | ormr2 | Multiplies a real matrix by the orthogonal matrix Q formed by gerq2. |
| :+1:  | ungr2 | Generates the complex unitary matrix Q of the RQ factorization formed by gerq2. |
| :+1:  | unmr2 | Multiplies a complex matrix by the unitary matrix Q formed by gerq2. |

#### Singular Value and Eigenvalue Problem Routines
| done?| name  | description             |
| ---- | ----- | ----------------------- |
| :+1: | gesvd | Computes the singular value decomposition of a general rectangular matrix.  |
| :+1: | heevd | Computes all eigenvalues and, optionally, all eigenvectors of a complex Hermitian matrix using divide and conquer algorithm. |
| :+1: | hegvd | Computes all eigenvalues and, optionally, all eigenvectors of a complex generalized Hermitian definite eigenproblem using divide and conquer algorithm. |
| :+1:  | heevr | Computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices. |
| f77  | heevx | Computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices. |
|      | gebrd | Reduces a general matrix to bidiagonal form.     |
|      | hetrd | Reduces a complex Hermitian matrix to tridiagonal form. |
|      | orgbr | Generates the real orthogonal matrix Q or PT determined by gebrd. |
|      | orgtr | Generates the real orthogonal matrix Q determined by sytrd. |
|      | ormtr | Multiplies a real matrix by the orthogonal matrix Q determined by sytrd. |
|      | syevd | Computes all eigenvalues and, optionally, all eigenvectors of a real symmetric matrix using divide and conquer algorithm. |
|      | sygvd | Computes all eigenvalues and, optionally, all eigenvectors of a real generalized symmetric definite eigenproblem using divide and conquer algorithm. |
|      | sytrd | Reduces a real symmetric matrix to tridiagonal form. |
|      | ungbr | Generates the complex unitary matrix Q or PT determined by gebrd. |
|      | ungtr | Generates the complex unitary matrix Q determined by hetrd. |
|      | unmtr | Multiplies a complex matrix by the unitary matrix Q determined by hetrd. |

##### Least squares
|done| name  | description                                    |
|----| ----- | ---------------------------------------------- |
|f77 | gels  | least squares using QR/LQ                      |
|f77 | gelst | least squares using QR/LQ with T matrix        |
|f77 | gelss | least squares using SVD, QR iteration          |
|f77 | gelsd | least squares using SVD, divide and conquer    |
|f77 | gelsy | least squares using complete orthogonal factor |
|f77 | getsls| least squares using tall-skinny QR/LQ          |
|f77 | gglse | equality-constrained least squares             |
|f77 | ggglm | Gauss-Markov linear model                      |

#### Other Auxiliary Routines

There are some other auxiliary lapack routines around, that may apear here:

| name      | Data Types | Description |
| --------- | ---------- | ------------|
| mfi_lartg | s, d, c, z | Generates a plane rotation with real cosine and real/complex sine. |

