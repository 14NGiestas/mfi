================================================
FILE: README.md
================================================
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
|  f77  | orgqr | Generates the real orthogonal matrix Q of the QR factorization formed by geqrf. |
|  f77  | orgrq | Generates the real orthogonal matrix Q of the RQ factorization formed by gerqf. |
|  f77  | ormqr | Multiplies a real matrix by the orthogonal matrix Q of the QR factorization formed by geqrf. |
|  f77  | ormrq | Multiplies a real matrix by the orthogonal matrix Q of the RQ factorization formed by gerqf. |
|       | sytrf | Computes the Bunch-Kaufman factorization of a symmetric matrix.                        |
|       | trtrs | Solves a system of linear equations with a triangular coefficient matrix, with multiple right-hand sides. |
|  f77  | ungqr | Generates the complex unitary matrix Q of the QR factorization formed by geqrf.  |
|  f77  | ungrq | Generates the complex unitary matrix Q of the RQ factorization formed by gerqf.  |
|  f77  | unmqr | Multiplies a complex matrix by the unitary matrix Q of the QR factorization formed by geqrf. |
|  f77  | unmrq | Multiplies a complex matrix by the unitary matrix Q of the RQ factorization formed by gerqf. |
|  f77  | org2r | |
|  f77  | orm2r | |
|  f77  | ung2r | |
|  f77  | unm2r | |
|  f77  | orgr2 | |
|  f77  | ormr2 | |
|  f77  | ungr2 | |
|  f77  | unmr2 | |

#### Singular Value and Eigenvalue Problem Routines
| done?| name  | description             |
| ---- | ----- | ----------------------- |
| :+1: | gesvd | Computes the singular value decomposition of a general rectangular matrix.  |
| :+1: | heevd | Computes all eigenvalues and, optionally, all eigenvectors of a complex Hermitian matrix using divide and conquer algorithm. |
| :+1: | hegvd | Computes all eigenvalues and, optionally, all eigenvectors of a complex generalized Hermitian definite eigenproblem using divide and conquer algorithm. |
| f77  | heevr | Computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices. |
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




================================================
FILE: common.fpp
================================================
#:mute

#:set REAL    = 'real(wp)'
#:set COMPLEX = 'complex(wp)'
#:set PREFIX = { &
    's': { 'type': 'real(wp)',    'wp': 'REAL32'}, &
    'd': { 'type': 'real(wp)',    'wp': 'REAL64'}, &
    'c': { 'type': 'complex(wp)', 'wp': 'REAL32'}, &
    'z': { 'type': 'complex(wp)', 'wp': 'REAL64'}, &
}

#:set ERROR = lambda pfx: { 'type': f'error: {pfx}', 'wp' : f'error: {pfx}' }

#:set mix    = lambda l, r: list(lp + rp for lp, rp in zip(l,r))
#:set split  = lambda pfx: list(pfx) if len(pfx) > 1 else pfx
#:set get_types = lambda pfxs: (pfxs[0], pfxs[0] if len(pfxs) == 1 else pfxs[1])
#:set get    = lambda pfx,what: PREFIX.get(pfx).get(what)
#:set prefix = lambda pfx, name: name.replace('?',pfx)
#:set kind   = lambda pfx: get(pfx,'wp')
#:set type   = lambda pfx: get(pfx,'type').replace('wp',kind(pfx))
#:set real   = lambda pfx: REAL.replace('wp',kind(pfx))
#:set complex= lambda pfx: COMPLEX.replace('wp',kind(pfx))

#:set SINGLE_TYPES  = ['s','c']
#:set DOUBLE_TYPES  = ['d','z']
#:set REAL_TYPES    = ['s','d']
#:set COMPLEX_TYPES = ['c','z']
#:set DEFAULT_TYPES = REAL_TYPES + COMPLEX_TYPES

#:set MIX_REAL_COMPLEX  = mix(REAL_TYPES,COMPLEX_TYPES)
#:set MIX_COMPLEX_REAL  = mix(COMPLEX_TYPES,REAL_TYPES)
#:set MIX_SINGLE_DOUBLE = mix(SINGLE_TYPES,DOUBLE_TYPES)
#:set MIX_DOUBLE_SINGLE = mix(DOUBLE_TYPES,SINGLE_TYPES)

#:def timeit(message, code)
block
real :: t1, t2
call cpu_time(t1)
$:code
call cpu_time(t2)
print '(A," (",G0,"s)")', ${message}$, t2-t1
end block
#:enddef

#:def random_number(type, name, shape='')
#:if type.startswith('complex')
    $:random_complex(type, name,shape)
#:else
    call random_number(${name}$)
#:endif
#:enddef

#:def random_complex(type, name, shape='')
#:set REAL = type.replace('complex','real')
block
    ${REAL}$ :: re${shape}$
    ${REAL}$ :: im${shape}$
    call random_number(im)
    call random_number(re)
    ${name}$ = cmplx(re,im)
end block
#:enddef

#! Handles parameters (usage: working precision)
#:def parameter(dtype, **kwargs)
#:for variable, value in kwargs.items()
    ${dtype}$, parameter :: ${variable}$ = ${value}$
#:endfor
#:enddef

#! Handles importing and setting precision constants in interfaces
#:def imports(pfxs)
#:set wps = set(list(map(kind, pfxs)))
#:if len(wps) > 1
    import :: ${', '.join(wps)}$
#:else
    import :: ${''.join(wps)}$
#:endif
#:enddef

#! Handles the input/output arguments
#:def args(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$) :: ${variable}$
#:endfor
#:enddef

#! Defines a optional variable, creating local corresponding variable by default
#:def optional(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$), optional :: ${variable}$
    ${dtype}$ :: local_${variable}$
#:endfor
#:enddef

#! Handles default values of a optional variable
#:def defaults(**kwargs)
#:for variable, default in kwargs.items()
    if (present(${variable}$)) then
        local_${variable}$ = ${variable}$
    else
        local_${variable}$ = ${default}$
    end if
#:endfor
#:enddef

#! Handles a value of "variable" depending on "condition"
#:def optval(condition, variable, true_value, false_value)
    if (${condition}$) then
        ${variable}$ = ${true_value}$
    else
        ${variable}$ = ${false_value}$
    end if
#:enddef

#:def interface(functions, procedure='procedure', name='')
interface ${name}$
    #:for function_name in functions
    ${procedure}$ :: ${function_name}$
    #:endfor
end interface
#:enddef

#! Interfaces for the original f77 routines
#! code must implement a routine interface
#:def f77_original(generic_name, prefixes, code)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set f90 = 'f77_' + prefix('',generic_name)
#:set f77 = [prefix(pfx,generic_name) for pfx in prefixes]
!> Generic old style interface for ${prefix('',generic_name).upper()}$.
!> Supports ${', '.join(prefixes)}$.
!> See also: [[${mfi}$]], ${'[[' + ']], [['.join(f77) + ']]'}$.
interface f77_${prefix('',generic_name)}$
#:for pfx in prefixes
#:set name = prefix(pfx,generic_name)
#:set pfxs = list(map(split,pfx))
!> Original interface for ${name.upper()}$
!> See also: [[${mfi}$]], [[${f90}$]].
$:code(name,pfxs)
#:endfor
end interface
#:enddef

#! Define a common interface with the original f77 interfaces
#! So you can call the original function without the prefix
#:def f77_improved(generic_name, prefixes)
#:set functions = map(lambda pfx: prefix(pfx,generic_name), prefixes)
$:interface(functions, name=f"f77_{prefix('',generic_name)}")
#:enddef

#! In case of missing functions / extensions you can pass a code
#! in which case must provide the routine implementation
#! Must be called inside a contains block
#:def f77_implement(generic_name, prefixes, code)
#:for pfx in prefixes
#:set name = prefix(pfx,generic_name)
#:set pfxs = list(map(split,pfx))
$:code(name,pfxs)
#:endfor
#:enddef

#:def mfi_interface(generic_name, prefixes)
#:set f77 = ['f77_' + prefix('',generic_name) + ':' + prefix(pfx,generic_name) for pfx in prefixes]
!> Generic modern interface for ${prefix('',generic_name).upper()}$.
!> Supports ${', '.join(prefixes)}$.
!> See also:
!> ${'[[' + ']], [['.join(f77) + ']]'}$.
#:set functions = map(lambda pfx: 'mfi_' + prefix(pfx,generic_name), prefixes)
$:interface(functions, &
            procedure='module procedure', &
            name=f"mfi_{prefix('',generic_name)}")
#:enddef

#! Implements the modern interface in code
#! for each supported prefix combination
#! Must be called inside a contains block
#:def mfi_implement(generic_name, prefixes, code)
#:for pfx in prefixes
#:set mfi_name  = 'mfi_' + prefix(pfx,generic_name)
#:set f77_name  =          prefix(pfx,generic_name)
#:set pfxs      = list(map(split,pfx))
#:set fun       = prefix('',generic_name)
!> Modern interface for [[f77_${fun}$:${f77_name}$]].
!> See also: [[mfi_${fun}$]], [[f77_${fun}$]].
$:code(mfi_name,f77_name,pfxs)
#:endfor
#:enddef


#! Implements the test for all interfaces
#! and each supported prefix combination
#! Must be called inside a contains block
#:def test_implement(generic_name, prefixes, code)
#:for pfx in prefixes
#:set f77  =          prefix(pfx,generic_name)
#:set f90 = 'f77_' + prefix('',generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfxs    = list(map(split,pfx))
$:code(f77,f90,mfi,pfxs)
#:endfor
#:enddef

#:def test_run(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
@:timeit("testing ${mfi}$ against ${f77}$", { call test_${f77}$ })
#:endfor
#:enddef

#:endmute



================================================
FILE: ford.md
================================================
---
project: MFI - Modern Fortran Interfaces 
src_dir: ./src
output_dir: ./api-reference
project_github: https://github.com/14NGiestas/mfi
project_website: 
summary: A collection of modern fortran interfaces for BLAS and LAPACK 
author: I. G. Pauli
author_description: 
github: https://github.org/14NGiestas
email: iangiestas@usp.br
macro: HAS_DECREMENT
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: true 
graph: true 
search: true 
sort: type
license: by-nc
max_frontpage_items: 4
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            iso_c_binding:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html#ISO_005fC_005fBINDING
---




================================================
FILE: fpm.toml
================================================
name    = "mfi"
version = "0.0.1"
license = "MIT"
author  = "Ian Giestas Pauli"
maintainer = "iangiestas@usp.br"

[library]
build-script="Makefile"

[build]
link = ["blas","lapack"]



================================================
FILE: LICENSE
================================================
MIT License

Copyright (c) 2020 Ian Giestas Pauli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



================================================
FILE: Makefile
================================================
FPP=fypp
#FYPPFLAGS=-DUBUNTU_WORKAROUND
fpp_files=$(shell find test src -name "*.fpp")
f90_files=$(patsubst %.fpp,%.f90,$(fpp_files))
all: $(f90_files)
%.f90: %.fpp; $(FPP) $(FYPPFLAGS) -I. $< $@
clean:; -rm $(f90_files)



================================================
FILE: src/f77/blas.fpp
================================================
#:mute
#:include "common.fpp"
#:include "src/f77/blas/asum_nrm2.fypp"
#:include "src/f77/blas/scal.fypp"
#:include "src/f77/blas/axpy.fypp"
#:include "src/f77/blas/copy_swap.fypp"
#:include "src/f77/blas/dot_product.fypp"
#:include "src/f77/blas/rot.fypp"
#:include "src/f77/blas/rotg.fypp"
#:include "src/f77/blas/rotm.fypp"
#:include "src/f77/blas/rotmg.fypp"
#:include "src/f77/blas/gbmv.fypp"
#:include "src/f77/blas/gemv.fypp"
#:include "src/f77/blas/ger_gerc_geru.fypp"
#:include "src/f77/blas/hbmv_sbmv.fypp"
#:include "src/f77/blas/hemv_symv.fypp"
#:include "src/f77/blas/her.fypp"
#:include "src/f77/blas/syr.fypp"
#:include "src/f77/blas/her_syr2.fypp"
#:include "src/f77/blas/hpmv_spmv.fypp"
#:include "src/f77/blas/hpr.fypp"
#:include "src/f77/blas/spr.fypp"
#:include "src/f77/blas/hpr_spr2.fypp"
#:include "src/f77/blas/tbmv_tbsv.fypp"
#:include "src/f77/blas/tpmv_tpsv.fypp"
#:include "src/f77/blas/trmv_trsv.fypp"
#:include "src/f77/blas/gemm.fypp"
#:include "src/f77/blas/hemm_symm.fypp"
#:include "src/f77/blas/herk.fypp"
#:include "src/f77/blas/syrk.fypp"
#:include "src/f77/blas/her2k.fypp"
#:include "src/f77/blas/syr2k.fypp"
#:include "src/f77/blas/trmm_trsm.fypp"
! BLAS Level 1 - Extensions
#:include "src/f77/blas/iamax_iamin.fypp"
#:include "src/f77/blas/iamin_stub.fypp"

#:set COLLECT = [                                              &
    ('?copy', DEFAULT_TYPES,                    copy_swap),    &
    ('?swap', DEFAULT_TYPES,                    copy_swap),    &
    ('?axpy', DEFAULT_TYPES,                    axpy),         &
    ('?dot',  REAL_TYPES,                       dot_product),  &
    ('?dotc', COMPLEX_TYPES,                    dot_product),  &
    ('?dotu', COMPLEX_TYPES,                    dot_product),  &
    ('?asum', REAL_TYPES + MIX_REAL_COMPLEX,    asum_nrm2),    &
    ('?nrm2', REAL_TYPES + MIX_REAL_COMPLEX,    asum_nrm2),    &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),          &
    ('?rotg', DEFAULT_TYPES,                    rotg),         &
    ('?rotm', REAL_TYPES,                       rotm),         &
    ('?rotmg',REAL_TYPES,                       rotmg),        &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),         &
    ('?gbmv', DEFAULT_TYPES,                    gbmv),         &
    ('?gemv', DEFAULT_TYPES,                    gemv),         &
    ('?ger',  REAL_TYPES,                       ger_gerc_geru),&
    ('?gerc', COMPLEX_TYPES,                    ger_gerc_geru),&
    ('?geru', COMPLEX_TYPES,                    ger_gerc_geru),&
    ('?hbmv', COMPLEX_TYPES,                    hbmv_sbmv),    &
    ('?hemv', COMPLEX_TYPES,                    hemv_symv),    &
    ('?her',  COMPLEX_TYPES,                    her),          &
    ('?her2', COMPLEX_TYPES,                    her_syr2),     &
    ('?hpmv', COMPLEX_TYPES,                    hpmv_spmv),    &
    ('?hpr',  COMPLEX_TYPES,                    hpr),          &
    ('?hpr2', COMPLEX_TYPES,                    hpr_spr2),     &
    ('?sbmv', REAL_TYPES,                       hbmv_sbmv),    &
    ('?spmv', REAL_TYPES,                       hpmv_spmv),    &
    ('?spr',  REAL_TYPES,                       spr),          &
    ('?spr2', REAL_TYPES,                       hpr_spr2),     &
    ('?symv', REAL_TYPES,                       hemv_symv),    &
    ('?syr',  REAL_TYPES,                       syr),          &
    ('?syr2', REAL_TYPES,                       her_syr2),     &
    ('?tbmv', DEFAULT_TYPES,                    tbmv_tbsv),    &
    ('?tbsv', DEFAULT_TYPES,                    tbmv_tbsv),    &
    ('?tpmv', DEFAULT_TYPES,                    tpmv_tpsv),    &
    ('?tpsv', DEFAULT_TYPES,                    tpmv_tpsv),    &
    ('?trmv', DEFAULT_TYPES,                    trmv_trsv),    &
    ('?trsv', DEFAULT_TYPES,                    trmv_trsv),    &
    ('?gemm', DEFAULT_TYPES,                    gemm),         &
    ('?hemm', COMPLEX_TYPES,                    hemm_symm),    &
    ('?herk', COMPLEX_TYPES,                    herk),         &
    ('?her2k',COMPLEX_TYPES,                    her2k),        &
    ('?symm', REAL_TYPES,                       hemm_symm),    &
    ('?syrk', REAL_TYPES,                       syrk),         &
    ('?syr2k',REAL_TYPES,                       syr2k),        &
    ('?trmm', DEFAULT_TYPES,                    trmm_trsm),    &
    ('?trsm', DEFAULT_TYPES,                    trmm_trsm),    &
]
#:endmute
!> Improved and original F77 interfaces for BLAS
module f77_blas
use iso_fortran_env
implicit none

#:for name, supported_types, code in COLLECT
$:f77_original(name, supported_types, code)
#:endfor

#:include "src/f77/blas/specific_interfaces.fypp"

! Extensions
! BLAS Level 1 - Utils / Extensions
#:if defined('MFI_EXTENSIONS')
  #:if defined('MFI_LINK_EXTERNAL')
! Link with a external source
$:f77_original('i?amax', DEFAULT_TYPES, iamax_iamin)
$:f77_original('i?amin', DEFAULT_TYPES, iamax_iamin)
  #:else
! Implement the blas extensions in
$:f77_improved('i?amax', DEFAULT_TYPES)
$:f77_improved('i?amin', DEFAULT_TYPES)
contains
$:f77_implement('i?amax', DEFAULT_TYPES, iamin_stub)
$:f77_implement('i?amin', DEFAULT_TYPES, iamin_stub)
  #:endif
#:endif

end module




================================================
FILE: src/f77/lapack.fpp
================================================
#:mute
#:include "common.fpp"
#:include "src/f77/lapack/lartg.fypp"
#:include "src/f77/lapack/geqrf_gerqf.fypp"
#:include "src/f77/lapack/getrf.fypp"
#:include "src/f77/lapack/getri.fypp"
#:include "src/f77/lapack/getrs.fypp"
#:include "src/f77/lapack/hetrf.fypp"
#:include "src/f77/lapack/gesvd.fypp"
#:include "src/f77/lapack/hegv.fypp"
#:include "src/f77/lapack/heevx.fypp"
#:include "src/f77/lapack/heevr.fypp"
#:include "src/f77/lapack/heevd.fypp"
#:include "src/f77/lapack/potrf_potri.fypp"
#:include "src/f77/lapack/potrs.fypp"
#:include "src/f77/lapack/pocon.fypp"
#:include "src/f77/lapack/gels_gelst_getsls.fypp"
#:include "src/f77/lapack/gelsd.fypp"
#:include "src/f77/lapack/gelss.fypp"
#:include "src/f77/lapack/gelsy.fypp"
#:include "src/f77/lapack/gglse.fypp"
#:include "src/f77/lapack/gglsm.fypp"
#:include "src/f77/lapack/org2r_orgr2_ung2r_ungr2.fypp"
#:include "src/f77/lapack/orgqr_orgrq_ungqr_ungrq.fypp"
#:include "src/f77/lapack/orm2r_ormr2_unm2r_unmr2.fypp"
#:include "src/f77/lapack/ormqr_ormrq_unmqr_unmrq.fypp"
#:set COLLECT = [                                  &
    ('?geqrf',  DEFAULT_TYPES, geqrf_gerqf),       &
    ('?gerqf',  DEFAULT_TYPES, geqrf_gerqf),       &
    ('?getrf',  DEFAULT_TYPES, getrf),             &
    ('?getri',  DEFAULT_TYPES, getri),             &
    ('?getrs',  DEFAULT_TYPES, getrs),             &
    ('?hetrf',  COMPLEX_TYPES, hetrf),             &
    ('?hegv',   COMPLEX_TYPES, hegv),              &
    ('?heevd',  COMPLEX_TYPES, heevd),             &
    ('?gesvd',  DEFAULT_TYPES, gesvd),             &
    ('?potrf',  DEFAULT_TYPES, potrf_potri),       &
    ('?potri',  DEFAULT_TYPES, potrf_potri),       &
    ('?potrs',  DEFAULT_TYPES, potrs),             &
    ('?pocon',  DEFAULT_TYPES, pocon),             &
    ('?heevx',  COMPLEX_TYPES, heevx),             &
    ('?heevr',  COMPLEX_TYPES, heevr),             &
    ('?gels',   DEFAULT_TYPES, gels_gelst_getsls), &
    ('?gelst',  DEFAULT_TYPES, gels_gelst_getsls), &
    ('?getsls', DEFAULT_TYPES, gels_gelst_getsls), &
    ('?gelsd',  DEFAULT_TYPES, gelsd),             &
    ('?gelss',  DEFAULT_TYPES, gelss),             &
    ('?gelsy',  DEFAULT_TYPES, gelsy),             &
    ('?gglse',  DEFAULT_TYPES, gglse),             &
    ('?gglsm',  DEFAULT_TYPES, gglsm),             &
    ('?org2r',  REAL_TYPES,    org2r_orgr2_ung2r_ungr2), &
    ('?orgr2',  REAL_TYPES,    org2r_orgr2_ung2r_ungr2), &
    ('?orm2r',  REAL_TYPES,    orm2r_ormr2_unm2r_unmr2), &
    ('?ormr2',  REAL_TYPES,    orm2r_ormr2_unm2r_unmr2), &
    ('?ormqr',  REAL_TYPES,    ormqr_ormrq_unmqr_unmrq), &
    ('?ormrq',  REAL_TYPES,    ormqr_ormrq_unmqr_unmrq), &
    ('?orgqr',  REAL_TYPES,    orgqr_orgrq_ungqr_ungrq), &
    ('?orgrq',  REAL_TYPES,    orgqr_orgrq_ungqr_ungrq), &
    ('?ung2r',  COMPLEX_TYPES, org2r_orgr2_ung2r_ungr2), &
    ('?ungr2',  COMPLEX_TYPES, org2r_orgr2_ung2r_ungr2), &
    ('?unm2r',  COMPLEX_TYPES, orm2r_ormr2_unm2r_unmr2), &
    ('?unmr2',  COMPLEX_TYPES, orm2r_ormr2_unm2r_unmr2), &
    ('?unmqr',  COMPLEX_TYPES, ormqr_ormrq_unmqr_unmrq), &
    ('?unmrq',  COMPLEX_TYPES, ormqr_ormrq_unmqr_unmrq), &
    ('?ungqr',  COMPLEX_TYPES, orgqr_orgrq_ungqr_ungrq), &
    ('?ungrq',  COMPLEX_TYPES, orgqr_orgrq_ungqr_ungrq), &
    ('?lartg',  DEFAULT_TYPES, lartg),             &
]
#:endmute
!> Improved and original F77 interfaces for LAPACK
module f77_lapack
use iso_fortran_env
implicit none

#:for name, supported_types, code in COLLECT
$:f77_original(name, supported_types, code)
#:endfor

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module




================================================
FILE: src/f77/blas/asum_nrm2.fypp
================================================
#:def asum_nrm2(NAME,pfxs)
#:set A, B = get_types(pfxs)
pure function ${NAME}$(n, x, incx)
$:imports(pfxs)
    ${real(A)}$ :: ${NAME}$
@:args(${type(B)}$, in, x(*))
@:args(integer,     in, n, incx)
end function
#:enddef



================================================
FILE: src/f77/blas/axpy.fypp
================================================
#:def axpy(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(n, a, x, incx, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(*), a)
@:args(${type(wp)}$, inout, y(*))
@:args(integer,  in,    n, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/copy_swap.fypp
================================================
#:def copy_swap(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(n, x, incx, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(*))
@:args(${type(wp)}$, inout, y(*))
@:args(integer,  in,    n, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/dot_product.fypp
================================================
#:def dot_product(NAME,pfxs)
#:set wp = pfxs[0]
pure function ${NAME}$(n, x, incx, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
    ${type(wp)}$ :: ${NAME}$
@:args(${type(wp)}$, in, x(*), y(*))
@:args(integer,  in, n, incx, incy)
end function
#:enddef



================================================
FILE: src/f77/blas/gbmv.fypp
================================================
#:def gbmv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*), x(*))
@:args(${type(wp)}$,  inout, y(*))
@:args(character, in,    trans)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, kl, ku, lda, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/gemm.fypp
================================================
#:def gemm(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*), b(ldb,*))
@:args(${type(wp)}$,  inout, c(ldc,*))
@:args(character, in,    transa, transb)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, k, lda, ldb, ldc)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/gemv.fypp
================================================
#:def gemv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*), x(*))
@:args(${type(wp)}$,  inout, y(*))
@:args(character, in,    trans)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, lda, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/ger_gerc_geru.fypp
================================================
#:def ger_gerc_geru(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(m, n, alpha, x, incx, y, incy, a, lda)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*), y(*))
@:args(${type(wp)}$,  inout, a(lda,*))
@:args(${type(wp)}$,  in,    alpha)
@:args(integer,   in,    m, n, lda, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/hbmv_sbmv.fypp
================================================
#:def hbmv_sbmv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*), x(*))
@:args(${type(wp)}$,  inout, y(*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/hemm_symm.fypp
================================================
#:def hemm_symm(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*), b(ldb,*))
@:args(${type(wp)}$,  inout, c(ldc,*))
@:args(character, in,    side, uplo)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, lda, ldb, ldc)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/hemv_symv.fypp
================================================
#:def hemv_symv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*), x(*))
@:args(${type(wp)}$,  inout, y(*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    n, lda, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/her.fypp
================================================
#:def her(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, a, lda)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*))
@:args(${type(wp)}$,  inout, a(lda,*))
@:args(character, in,    uplo)
@:args(real(wp),  in,    alpha)
@:args(integer,   in,    n, lda, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/her2k.fypp
================================================
#:def her2k(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  in,    b(ldb,*))
@:args(${type(wp)}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(${type(wp)}$,  in,    alpha)
@:args(real(wp),  in,    beta)
@:args(integer,   in,    n, k, lda, ldb, ldc)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/her_syr2.fypp
================================================
#:def her_syr2(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*), y(*))
@:args(${type(wp)}$,  inout, a(lda,*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha)
@:args(integer,   in,    n, lda, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/herk.fypp
================================================
#:def herk(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(real(wp),  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldc)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/hpmv_spmv.fypp
================================================
#:def hpmv_spmv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    ap(*), x(*))
@:args(${type(wp)}$,  inout, y(*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    n, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/hpr.fypp
================================================
#:def hpr(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, ap)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*))
@:args(${type(wp)}$,  inout, ap(*))
@:args(character, in,    uplo)
@:args(real(wp),  in,    alpha)
@:args(integer,   in,    n, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/hpr_spr2.fypp
================================================
#:def hpr_spr2(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, y, incy, ap)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*), y(*))
@:args(${type(wp)}$,  inout, ap(*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha)
@:args(integer,   in,    n, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/iamax_iamin.fypp
================================================
#:def iamax_iamin(NAME,pfxs)
#:set wp = pfxs[0]
pure function ${NAME}$(n, x, incx)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
    integer :: ${NAME}$
@:args(${type(wp)}$, in, x(*))
@:args(integer,  in, n, incx)
end function
#:enddef



================================================
FILE: src/f77/blas/iamin_stub.fypp
================================================
#:def iamin_stub(NAME,pfxs)
#:set wp = pfxs[0]
pure function ${NAME}$(n, x, incx)
@:parameter(integer, wp=${kind(wp)}$)
    integer :: ${NAME}$
@:args(${type(wp)}$, in, x(*))
@:args(integer,  in, n, incx)
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        ${NAME}$ = 0
        return
    end if
#:if type(wp) == complex(wp)
    ${NAME}$ = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
#:else
    ${NAME}$ = minloc(x(1:n:incx),dim=1)
#:endif
end function
#:enddef



================================================
FILE: src/f77/blas/rot.fypp
================================================
#:mute
subroutine {s,d,c,z}rot (
    integer  n,
    type(wp) x(*),
    integer  incx,
    type(wp) y(*),
    integer  incy,
    real(wp) c,
    type(wp) s
)
subroutine {sc,zd}rot (
    integer  n,
    type(wp) x(*),
    integer  incx,
    type(wp) y(*),
    integer  incy,
    real(wp) c,
    real(wp) s
)
#:def rot(name,pfxs)
#:set A, B = get_types(pfxs)
!> ${name.upper()}$ applies a plane rotation.
pure subroutine ${name}$(n, x, incx, y, incy, c, s)
$:imports(pfxs)
@:args(${type(A)}$,   in, x(*), y(*))
@:args(integer,       in, n, incx, incy)
@:args(${real(A)}$,   in, c)
@:args(${type(B)}$,   in, s)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/blas/rotg.fypp
================================================
#:mute
subroutine {s,d,c,z}rotg (
    type(wp)  a,
    type(wp)  b,
    real(wp)  c,
    type(wp)  s
)
#:def rotg(NAME,pfxs)
#:set wp = pfxs[0]
!>${NAME}$ generates a Givens rotation with real cosine and complex sine:
#:if type(wp) == real(wp)
!>```
!> [  c  s ] [ a ] = [ r ]
!> [ -s  c ] [ b ]   [ 0 ]
!>```
!> satisfying `c**2 + s**2 = 1`.
#:elif type(wp) == complex(wp)
!>```
!>  [  c         s ] [ a ] = [ r ]
!>  [ -conjg(s)  c ] [ b ]   [ 0 ]
!>```
!> where c is real, s is complex, and `c**2 + conjg(s)*s = 1`.
#:endif
pure subroutine ${NAME}$(a, b, c, s)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a, b)
@:args(${real(wp)}$, out, c)
@:args(${type(wp)}$,      out, s)
end subroutine

#:enddef
#:endmute



================================================
FILE: src/f77/blas/rotm.fypp
================================================
#:def rotm(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(n, x, incx, y, incy, param)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, x(*), y(*))
@:args(${type(wp)}$, in, param(5))
@:args(integer,  in, n, incx, incy)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/rotmg.fypp
================================================
#:def rotmg(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(d1, d2, x1, y1, param)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    y1)
@:args(${type(wp)}$, out,   param(5))
@:args(${type(wp)}$, inout, d1, d2, x1)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/scal.fypp
================================================
#:def scal(NAME,pfxs)
#:set A, B = get_types(pfxs)
!> ${NAME.upper()}$ scales a vector by a constant.
pure subroutine ${NAME}$(n, a, x, incx)
$:imports(pfxs)
@:args(${type(A)}$, inout, x(*))
@:args(${type(B)}$, in,    a)
@:args(integer,  in,    n, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/specific_interfaces.fypp
================================================
!> ?lamch supports s, d. See [[mfi_lamch]] for the modern version.
interface
    !> SLAMCH determines single precision machine parameters.
    pure real(REAL32) function slamch(cmach)
        import :: REAL32
        character, intent(in) :: cmach
    end function

    !> DLAMCH determines double precision machine parameters.
    pure real(REAL64) function dlamch(cmach)
        import :: REAL64
        character, intent(in) :: cmach
    end function
end interface

interface
    !> Compute the inner product of two vectors with extended
    !> precision accumulation.
    !>
    !> Returns S.P. result with dot product accumulated in D.P.
    !> SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
    !> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
    !> defined in a similar way using INCY.
    pure function sdsdot(n, sb, sx, incx, sy, incy)
        import :: REAL32
        integer, parameter :: wp = REAL32
        real(wp) :: sdsdot
        real(wp), intent(in) :: sx(*)
        real(wp), intent(in) :: sy(*)
        real(wp), intent(in) :: sb
        integer, intent(in) :: n
        integer, intent(in) :: incx
        integer, intent(in) :: incy
    end function

    !> Compute the inner product of two vectors with extended
    !> precision accumulation and result.
    !>
    !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
    !> DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
    !> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
    !> defined in a similar way using INCY.
    pure function dsdot(n, sx, incx, sy, incy)
        import :: REAL32, REAL64
        integer, parameter :: sp = REAL32
        integer, parameter :: dp = REAL64
        real(dp) :: dsdot
        real(sp), intent(in) :: sx(*)
        real(sp), intent(in) :: sy(*)
        integer,  intent(in) :: n
        integer,  intent(in) :: incx
        integer,  intent(in) :: incy
    end function
end interface



================================================
FILE: src/f77/blas/spr.fypp
================================================
#:def spr(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, ap)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*))
@:args(${type(wp)}$,  inout, ap(*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha)
@:args(integer,   in,    n, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/syr.fypp
================================================
#:def syr(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, a, lda)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    x(*))
@:args(${type(wp)}$,  inout, a(lda,*))
@:args(character, in,    uplo)
@:args(${type(wp)}$,  in,    alpha)
@:args(integer,   in,    n, lda, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/syr2k.fypp
================================================
#:def syr2k(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  in,    b(ldb,*))
@:args(${type(wp)}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldb, ldc)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/syrk.fypp
================================================
#:def syrk(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(${type(wp)}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldc)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/tbmv_tbsv.fypp
================================================
#:def tbmv_tbsv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  inout, x(*))
@:args(character, in,    uplo, trans, diag)
@:args(integer,   in,    n, k, lda, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/tpmv_tpsv.fypp
================================================
#:def tpmv_tpsv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, diag, n, ap, x, incx)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    ap(*))
@:args(${type(wp)}$,  inout, x(*))
@:args(character, in,    uplo, trans, diag)
@:args(integer,   in,    n, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/trmm_trsm.fypp
================================================
#:def trmm_trsm(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  inout, b(ldb,*))
@:args(character, in,    side, uplo, transa, diag)
@:args(${type(wp)}$,  in,    alpha)
@:args(integer,   in,    m, n, lda, ldb)
end subroutine
#:enddef



================================================
FILE: src/f77/blas/trmv_trsv.fypp
================================================
#:def trmv_trsv(NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${NAME}$(uplo, trans, diag, n, a, lda, x, incx)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in,    a(lda,*))
@:args(${type(wp)}$,  inout, x(*))
@:args(character, in,    uplo, trans, diag)
@:args(integer,   in,    n, lda, incx)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/gels_gelst_getsls.fypp
================================================
#:mute
subroutine {s,d,c,z}gels / gelst / getsls (
        character   trans,
        integer         m,
        integer         n,
        integer      nrhs,
        type(wp) a(lda,*),
        integer       lda,
        type(wp) b(ldb,*),
        integer       ldb,
        type(wp)  work(*),
        integer     lwork,
        integer      info
)
#:def gels_gelst_getsls(NAME,pfxs)
#:set wp=pfxs[0]
!> ${NAME.upper()}$ solves overdetermined or underdetermined systems for GE matrices
#:if NAME.endswith('t')
!> using QR or LQ factorization with compact WY representation of Q.
#:endif
pure subroutine ${NAME}$(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(character,   in, trans)
@:args(${type(wp)}$, inout, a(lda,*), b(ldb,*))
@:args(${type(wp)}$,   out, work(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, lda, ldb, nrhs, lwork)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/gelsd.fypp
================================================
#:mute

subroutine {s,d,c,z}gelsd(
        integer  in             m,
        integer  in             n,
        integer  in          nrhs,
        type(wp) inout  a(lda,*),
        integer  in           lda,
        type(wp) inout  b(ldb,*),
        integer  in           ldb,
        type(wp) out         s(*),
        type(wp) in         rcond,
        integer  out         rank,
        type(wp) out      work(*),
        integer  in         lwork,
        integer  out     iwork(*),
        integer  out         info
)
#:def gelsd(NAME,pfxs)
#:set wp=pfxs[0]
!> ${NAME.upper()}$ computes the minimum-norm solution to a linear least squares problem for GE matrices
pure subroutine ${NAME}$(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    rcond)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, inout, b(ldb,*))
@:args(${type(wp)}$, out,   s(*))
@:args(${type(wp)}$, out,   work(*))
@:args(integer, out, iwork(*))
@:args(integer, out, info, rank)
@:args(integer, in, n, m, nrhs, lda, ldb, lwork)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/gelss.fypp
================================================
#:mute

subroutine {s,d,c,z}gelss(
        integer  in             m,
        integer  in             n,
        integer  in          nrhs,
        type(wp) inout  a(lda,*),
        integer  in           lda,
        type(wp) inout  b(ldb,*),
        integer  in           ldb,
        type(wp) out         s(*),
        type(wp) in         rcond,
        integer  out         rank,
        type(wp) out      work(*),
        integer  in         lwork,
        integer  out         info
)
#:def gelss(NAME,pfxs)
#:set wp=pfxs[0]
!> ${NAME.upper()}$ solves overdetermined or underdetermined systems for GE matrices
pure subroutine ${NAME}$(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    rcond)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, inout, b(ldb,*))
@:args(${type(wp)}$, out,   s(*))
@:args(${type(wp)}$, out,   work(*))
@:args(integer, out, info, rank)
@:args(integer, in, n, m, nrhs, lda, ldb, lwork)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/gelsy.fypp
================================================
#:mute

subroutine {s,d,c,z}gelsy(
        integer  in             m,
        integer  in             n,
        integer  in          nrhs,
        type(wp) inout  a(lda,*),
        integer  in           lda,
        type(wp) inout  b(ldb,*),
        integer  in           ldb,
        integer  inout    jpvt(*),
        type(wp) in         rcond,
        integer  out         rank,
        type(wp) out      work(*),
        integer  in         lwork,
        integer  out         info
)
#:def gelsy(NAME,pfxs)
#:set wp=pfxs[0]
!> ${NAME.upper()}$ solves overdetermined or underdetermined systems for GE matrices
pure subroutine ${NAME}$(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    rcond)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, inout, b(ldb,*))
@:args(integer,  inout,  jpvt(*))
@:args(${type(wp)}$, out,   work(*))
@:args(integer, out, info, rank)
@:args(integer, in, n, m, nrhs, lda, ldb, lwork)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/geqrf_gerqf.fypp
================================================
#:def geqrf_gerqf(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(m,n,a,lda,tau,work,lwork,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$,   out, tau(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, lda, lwork)
@:args(${type(wp)}$, inout, work(*))
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/gesv.fypp
================================================
#:mute
subroutine {s,d,c,z}gesv (
    integer  in           n,
    integer  in        nrhs,
    type(wp) inout a(lda,*),
    integer  in         lda,
    integer  out    ipiv(*),
    type(wp) inout b(ldb,*),
    integer  in         ldb,
    integer  out       info
)

#:def gesv(NAME,pfxs)
#:set wp=pfxs[0]
!> ${NAME.upper()}$ computes the solution to system
!> of linear equations \( A \times X = B \) for GE matrices
pure subroutine ${NAME}$(n,nhrs,a,lda,ipiv,b,ldb,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  inout, a(lda,*), b(ldb,*))
@:args(integer,   out,   ipiv(*))
@:args(integer,   out,   info)
@:args(integer,   in,    n, nrhs, lda, ldb)
end subroutine
#:enddef

subroutine dsgesv  (
    integer   n,
    integer   nrhs,
    real(dp)  a( lda, * ),
    integer   lda,
    integer   ipiv( * ),
    real(dp)  b( ldb, * ),
    integer   ldb,
    real(dp)  x( ldx, * ),
    integer   ldx,
    real(dp)  work( n, * ),
    real(sp)  swork( * ),
    integer   iter,
    integer   info
)

#:def gesv_mixed(NAME,pfxs)
#:set wp=pfxs[0]
!> ${NAME.upper()}$ computes the solution to system of linear equations \( A \times X = B \) for GE matrices (mixed precision with iterative refinement)
pure subroutine ${NAME}$(n,nhrs,a,lda,ipiv,b,ldb,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  inout, a(lda,*), b(ldb,*))
@:args(integer,   out,   ipiv(*))
@:args(integer,   out,   info)
@:args(integer,   in,    n, nrhs, lda, ldb)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/gesvd.fypp
================================================
#:def gesvd(NAME,pfxs)
#:set wp=pfxs[0]
#:if type(wp) == complex(wp)
pure subroutine ${NAME}$(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
#:else
pure subroutine ${NAME}$(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
#:endif
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  inout, a(lda,*))
@:args(${real(wp)}$, out, s(*))
@:args(${type(wp)}$,  out,   u(ldu,*), vt(ldvt,*))
@:args(integer,   out,   info)
@:args(character, in,    jobu, jobvt)
@:args(integer,   in,    m, n, lda, ldu, ldvt, lwork)
@:args(${type(wp)}$,  inout, work(*))
#:if type(wp) == complex(wp)
@:args(${real(wp)}$, in, rwork(*))
#:endif
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/getrf.fypp
================================================
#:def getrf(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(m,n,a,lda,ipiv,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(integer,    out, ipiv(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, lda)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/getri.fypp
================================================
#:def getri(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(n,a,lda,ipiv,work,lwork,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, inout, work(*))
@:args(integer,     in, ipiv(*))
@:args(integer,    out, info)
@:args(integer,     in, n, lda, lwork)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/getrs.fypp
================================================
#:def getrs(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, inout, b(ldb,*))
@:args(character,   in, trans)
@:args(integer,     in, ipiv(*))
@:args(integer,    out, info)
@:args(integer,     in, n, nrhs, lda, ldb)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/gglse.fypp
================================================
#:mute
subroutine {s,d,c,z}gglse (
    integer  in    m
    integer  in    n
    integer  in    p
    type(wp) inout a(lda,*)
    integer  in    lda
    type(wp) inout b(ldb,*)
    integer  in    ldb
    type(wp) inout c(*)
    type(wp) inout d(*)
    type(wp) out   x(*)
    type(wp) out   work(*)
    integer  in    lwork
    integer  out   info
)
#:def gglse(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*), b(ldb,*), c(*), d(*))
@:args(${type(wp)}$,   out, work(*), x(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, p, lda, ldb, lwork)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/gglsm.fypp
================================================
#:mute
subroutine {s,d,c,z}gglsm (
    integer  in    n
    integer  in    m
    integer  in    p
    type(wp) inout a(lda,*)
    integer  in    lda
    type(wp) inout b(ldb,*)
    integer  in    ldb
    type(wp) inout d(*)
    type(wp) out   x(*)
    type(wp) out   y(*)
    type(wp) out   work(*)
    integer  in    lwork
    integer  out   info
)

#:def gglsm(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*), b(ldb,*), d(*))
@:args(${type(wp)}$,   out, work(*), x(*), y(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, p, lda, ldb, lwork)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/heevd.fypp
================================================
#:def heevd(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(lda,*))
@:args(${real(wp)}$, out,   w(*))
@:args(integer,       out,   info)
@:args(character,     in,    jobz, uplo)
@:args(integer,       in,    n, lda, lwork, lrwork, liwork)
@:args(${type(wp)}$,      inout, work(*))
@:args(${real(wp)}$, inout, rwork(*))
@:args(integer,       inout, iwork(*))
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/heevr.fypp
================================================
#:def heevr(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,&
                         isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(lda,*), z(ldz, *))
@:args(${real(wp)}$, out,   w(*))
@:args(integer,       out,   info)
@:args(character,     in,    jobz, uplo, range)
@:args(${real(wp)}$, in,    vl, vu, abstol)
@:args(integer,       in,    n, m, lda, ldz, il, iu, lwork, lrwork, liwork)
@:args(integer,       in,    isuppz(*))
@:args(${type(wp)}$,      inout, work(*))
@:args(${real(wp)}$, inout, rwork(*))
@:args(integer,       inout, iwork(*))
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/heevx.fypp
================================================
#:def heevx(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,&
                         work,lwork,rwork,lrwork,iwork,liwork,ifail,info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(lda,*), z(ldz, *))
@:args(${real(wp)}$, out,   w(*))
@:args(integer,       out,   info)
@:args(character,     in,    jobz, uplo, range)
@:args(${real(wp)}$, in,    vl, vu, abstol)
@:args(integer,       in,    n, m, lda, ldz, il, iu, lwork, lrwork, liwork, ifail)
@:args(${type(wp)}$,      inout, work(*))
@:args(${real(wp)}$, inout, rwork(*))
@:args(integer,       inout, iwork(*))
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/hegv.fypp
================================================
#:def hegv(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  inout, a(lda,*))
@:args(${type(wp)}$,  inout, b(ldb,*))
@:args(${real(wp)}$, out, w(*))
@:args(integer,   out,   info)
@:args(character, in,    jobz, uplo)
@:args(integer,   in,    n, itype, lda, ldb, lwork)
@:args(${type(wp)}$,  inout, work(*))
@:args(${real(wp)}$, in, rwork(*))
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/hetrf.fypp
================================================
#:def hetrf(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(uplo, n, a, lda, ipiv, work, lwork, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(character,   in, uplo)
@:args(integer,     in, ipiv(*))
@:args(${type(wp)}$, inout, work(*))
@:args(integer,    out, info)
@:args(integer,     in, n, lda, lwork)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/lartg.fypp
================================================
#:def lartg(NAME,pfxs)
#:set wp=pfxs[0]
#:set wp = pfxs[0]
pure subroutine ${NAME}$(f, g, c, s, r)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(real(wp), inout, c)
@:args(${type(wp)}$, inout, f, g, r, s)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/org2r_orgr2_ung2r_ungr2.fypp
================================================
#:mute
subroutine {s,d}org{2r,r2}(
subroutine {c,z}ung{2r,r2}(
        integer m,
        integer n,
        integer k,
        type(wp) a(lda,*),
        integer lda,
        type(wp) tau(*),
        type(wp) work(*),
        integer info
)

#:def org2r_orgr2_ung2r_ungr2(NAME,pfxs)
#:set wp=pfxs[0]
#:if NAME.endswith('2r')
!> This routine generates an \(M \times N \) ${type(wp)}$
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:${wp}$geqrf]].
#:elif NAME.endswith('r2')
!> This routine generates an \(M \times N \) ${type(wp)}$
!> matrix \( Q \) with orthonormal rows,
!> which is defined as the last \( M \) rows of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:${wp}$gerqf]].
#:endif
pure subroutine ${NAME}$(m, n, k, a, lda, tau, work, info)
$:imports(pfxs)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, out,   work(*))
@:args(${type(wp)}$, in,    tau(*))
@:args(integer,      in,    m, n, k, lda)
@:args(integer,      out,   info)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/orgqr_orgrq_ungqr_ungrq.fypp
================================================
#:mute
subroutine {s,d}org{qr,rq}(
subroutine {c,z}ung{qr,rq}(
        integer m,
        integer n,
        integer k,
        type(wp) a(lda,*),
        integer lda,
        type(wp) tau(*),
        type(wp) work(*),
        integer lwork,
        integer info
)

#:def orgqr_orgrq_ungqr_ungrq(NAME,pfxs)
#:set wp=pfxs[0]
!> This routine generates an \(M \times N \) ${type(wp)}$
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
#:if NAME.endswith('qr')
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:${wp}$geqrf]].
#:elif NAME.endswith('rq')
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger . . . H(k)^\dagger \)
!> as returned by [[f77_gerqf:${wp}$gerqf]].
#:endif
pure subroutine ${NAME}$(m, n, k, a, lda, tau, work, lwork, info)
$:imports(pfxs)
@:args(${type(wp)}$, inout, a(lda,*))
@:args(${type(wp)}$, out,   work(*))
@:args(${type(wp)}$, in,    tau(*))
@:args(integer,      in,    m, n, k, lda, lwork)
@:args(integer,      out,   info)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/orm2r_ormr2_unm2r_unmr2.fypp
================================================
#:mute
subroutine {s,d}orm{2r,r2} (
subroutine {c,z}unm{2r,r2} (
    character side,
    character trans,
    integer  m,
    integer  n,
    integer  k,
    type(wp) a(lda,*),
    integer  lda,
    type(wp) tau(*),
    type(wp) c(ldc,*),
    integer  ldc,
    type(wp) work(*),
    integer  info
)

#:def orm2r_ormr2_unm2r_unmr2(NAME,pfxs)
#:set wp=pfxs[0]
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
#:if NAME.endswith('2r')
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:${wp}$geqrf]].
#:elif NAME.endswith('r2')
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:${wp}$gerqf]].
#:endif
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine ${NAME}$(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
$:imports(pfxs)
@:args(character,    in,    side, trans)
@:args(${type(wp)}$, inout, a(lda,*), c(ldc,*))
@:args(${type(wp)}$, out,   work(*))
@:args(${type(wp)}$, in,    tau(*))
@:args(integer,      in,    m, n, k, lda, ldc)
@:args(integer,      out,   info)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/ormqr_ormrq_unmqr_unmrq.fypp
================================================
#:mute
subroutine {s,d}orm{qr,rq} (
subroutine {c,z}unm{qr,rq} (
    character side,
    character trans,
    integer  m,
    integer  n,
    integer  k,
    type(wp) a(lda,*),
    integer  lda,
    type(wp) tau(*),
    type(wp) c(ldc,*),
    integer  ldc,
    type(wp) work(*),
    integer  lwork,
    integer  info
)

#:def ormqr_ormrq_unmqr_unmrq(NAME,pfxs)
#:set wp=pfxs[0]
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
#:if NAME.endswith('qr')
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:${wp}$geqrf]].
#:elif NAME.endswith('rq')
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:${wp}$gerqf]].
#:endif
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine ${NAME}$(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
$:imports(pfxs)
@:args(character,    in,    side, trans)
@:args(${type(wp)}$, inout, a(lda,*), c(ldc,*))
@:args(${type(wp)}$, out,   work(*))
@:args(${type(wp)}$, in,    tau(*))
@:args(integer,      in,    m, n, k, lda, ldc, lwork)
@:args(integer,      out,   info)
end subroutine
#:enddef
#:endmute



================================================
FILE: src/f77/lapack/pocon.fypp
================================================
#:mute

!subroutine ?pocon (
!   character                      uplo,
!   integer                           n,
!   type(wp), dimension( lda, * )     a,
!   integer                         lda,
!   real(wp)                      anorm,
!   real(wp)                      rcond,
!   type(wp),    dimension( * )    work,
! if real(wp)S then,
!   integer,      dimension( * )   iwork,
! if complex(wp)S then,
!   real(wp),     dimension( * )   rwork,
!   integer                        info
!)
#:def pocon(NAME,pfxs)
#:set wp=pfxs[0]
#:set wp = pfxs[0]
!> ${NAME}$ estimates the reciprocal of the condition number (in the
!> 1-norm) of a ${type(wp)}$ Hermitian positive definite matrix using the
!> Cholesky factorization \( A = U^\dagger U \) or \( A = LL^\dagger |) computed by ${wp.upper()}$POTRF.
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
#:if type(wp) == real(wp)
pure subroutine ${NAME}$(uplo, n, a, lda, anorm, rcond, work, iwork, info)
#:elif type(wp) == complex(wp)
pure subroutine ${NAME}$(uplo, n, a, lda, anorm, rcond, work, rwork, info)
#:endif
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(character,     in,    uplo)
@:args(integer,       in,    n, lda)
@:args(${type(wp)}$,      inout, a(lda,*))
@:args(${real(wp)}$, in,    anorm)
@:args(${real(wp)}$, out,   rcond)
@:args(${type(wp)}$,      inout, work(*))
#:if type(wp) == real(wp)
@:args(integer,       inout, iwork(*))
#:elif type(wp) == complex(wp)
@:args(${real(wp)}$, inout, rwork(*))
#:endif
@:args(integer,       out,   info)
end subroutine
#:enddef

#:endmute



================================================
FILE: src/f77/lapack/potrf_potri.fypp
================================================
#:def potrf_potri(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(uplo, n, a, lda, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in, a(lda,*))
@:args(character, in,  uplo)
@:args(integer,   in,  n, lda)
@:args(integer,   out, info)
end subroutine
#:enddef



================================================
FILE: src/f77/lapack/potrs.fypp
================================================
#:def potrs(NAME,pfxs)
#:set wp=pfxs[0]
pure subroutine ${NAME}$(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: ${kind(wp)}$
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  in, a(lda,*))
@:args(${type(wp)}$,  in, b(ldb,*))
@:args(character, in,  uplo)
@:args(integer,   in,  n, nrhs, lda, ldb)
@:args(integer,   out, info)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas.fpp
================================================
#:mute
#:include "common.fpp"
#:include "src/mfi/blas/lamch.fypp"
#:include "src/mfi/blas/asum_nrm2.fypp"
#:include "src/mfi/blas/axpy.fypp"
#:include "src/mfi/blas/scal.fypp"
#:include "src/mfi/blas/copy_swap.fypp"
#:include "src/mfi/blas/dot_product.fypp"
#:include "src/mfi/blas/rot.fypp"
#:include "src/mfi/blas/rotm.fypp"
#:include "src/mfi/blas/iamin_iamax.fypp"
#:include "src/mfi/blas/gbmv.fypp"
#:include "src/mfi/blas/gemv.fypp"
#:include "src/mfi/blas/ger_gerc_geru.fypp"
#:include "src/mfi/blas/hbmv_sbmv.fypp"
#:include "src/mfi/blas/hemv_symv.fypp"
#:include "src/mfi/blas/her.fypp"
#:include "src/mfi/blas/syr.fypp"
#:include "src/mfi/blas/her_syr2.fypp"
#:include "src/mfi/blas/hpmv_spmv.fypp"
#:include "src/mfi/blas/hpr.fypp"
#:include "src/mfi/blas/spr.fypp"
#:include "src/mfi/blas/hpr_spr2.fypp"
#:include "src/mfi/blas/tbmv_tbsv.fypp"
#:include "src/mfi/blas/tpmv_tpsv.fypp"
#:include "src/mfi/blas/trmv_trsv.fypp"
#:include "src/mfi/blas/gemm.fypp"
#:include "src/mfi/blas/hemm_symm.fypp"
#:include "src/mfi/blas/herk.fypp"
#:include "src/mfi/blas/syrk.fypp"
#:include "src/mfi/blas/her2k.fypp"
#:include "src/mfi/blas/syr2k.fypp"
#:include "src/mfi/blas/trmm_trsm.fypp"
#:set COLLECT = [                                              &
    ('?copy', DEFAULT_TYPES,                    copy_swap),    &
    ('?swap', DEFAULT_TYPES,                    copy_swap),    &
    ('?axpy', DEFAULT_TYPES,                    axpy),         &
    ('?dot',  REAL_TYPES,                       dot_product),  &
    ('?dotc', COMPLEX_TYPES,                    dot_product),  &
    ('?dotu', COMPLEX_TYPES,                    dot_product),  &
    ('?asum', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),    &
    ('?nrm2', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),    &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),          &
    ('?rotm', REAL_TYPES,                       rotm),         &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),         &
    ('?gbmv', DEFAULT_TYPES,                    gbmv),         &
    ('?gemv', DEFAULT_TYPES,                    gemv),         &
    ('?ger',  REAL_TYPES,                       ger_gerc_geru),&
    ('?gerc', COMPLEX_TYPES,                    ger_gerc_geru),&
    ('?geru', COMPLEX_TYPES,                    ger_gerc_geru),&
    ('?hbmv', COMPLEX_TYPES,                    hbmv_sbmv),    &
    ('?hemv', COMPLEX_TYPES,                    hemv_symv),    &
    ('?her',  COMPLEX_TYPES,                    her),          &
    ('?her2', COMPLEX_TYPES,                    her_syr2),     &
    ('?hpmv', COMPLEX_TYPES,                    hpmv_spmv),    &
    ('?hpr',  COMPLEX_TYPES,                    hpr),          &
    ('?hpr2', COMPLEX_TYPES,                    hpr_spr2),     &
    ('?sbmv', REAL_TYPES,                       hbmv_sbmv),    &
    ('?spmv', REAL_TYPES,                       hpmv_spmv),    &
    ('?spr',  REAL_TYPES,                       spr),          &
    ('?spr2', REAL_TYPES,                       hpr_spr2),     &
    ('?symv', REAL_TYPES,                       hemv_symv),    &
    ('?syr',  REAL_TYPES,                       syr),          &
    ('?syr2', REAL_TYPES,                       her_syr2),     &
    ('?tbmv', DEFAULT_TYPES,                    tbmv_tbsv),    &
    ('?tbsv', DEFAULT_TYPES,                    tbmv_tbsv),    &
    ('?tpmv', DEFAULT_TYPES,                    tpmv_tpsv),    &
    ('?tpsv', DEFAULT_TYPES,                    tpmv_tpsv),    &
    ('?trmv', DEFAULT_TYPES,                    trmv_trsv),    &
    ('?trsv', DEFAULT_TYPES,                    trmv_trsv),    &
    ('?gemm', DEFAULT_TYPES,                    gemm),         &
    ('?hemm', COMPLEX_TYPES,                    hemm_symm),    &
    ('?herk', COMPLEX_TYPES,                    herk),         &
    ('?her2k',COMPLEX_TYPES,                    her2k),        &
    ('?symm', REAL_TYPES,                       hemm_symm),    &
    ('?syrk', REAL_TYPES,                       syrk),         &
    ('?syr2k',REAL_TYPES,                       syr2k),        &
    ('?trmm', DEFAULT_TYPES,                    trmm_trsm),    &
    ('?trsm', DEFAULT_TYPES,                    trmm_trsm),    &
    ('?lamch',REAL_TYPES,                       lamch),        &
]
#:endmute
!> Modern fortran interfaces for BLAS
module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotg  => f77_rotg
use f77_blas, only: mfi_rotmg => f77_rotmg
implicit none

#:for name, supported_types, code in COLLECT
$:mfi_interface(name, supported_types)
#:endfor

! Extensions
! BLAS level 1 - Utils / Extensions
#:if defined('MFI_EXTENSIONS')
$:mfi_interface('i?amax', DEFAULT_TYPES)
$:mfi_interface('i?amin', DEFAULT_TYPES)
#:endif

contains


#:for name, supported_types, code in COLLECT
$:mfi_implement(name, supported_types, code)
#:endfor

! Extensions
! BLAS level 1 - Utils / Extensions
#:if defined('MFI_EXTENSIONS')
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)
#:endif

end module



================================================
FILE: src/mfi/lapack.fpp
================================================
#:mute
#:include "common.fpp"
#:include "src/mfi/lapack/geqrf_gerqf.fypp"
#:include "src/mfi/lapack/getrf.fypp"
#:include "src/mfi/lapack/getri.fypp"
#:include "src/mfi/lapack/getrs.fypp"
#:include "src/mfi/lapack/hetrf.fypp"
#:include "src/mfi/lapack/gesvd.fypp"
#:include "src/mfi/lapack/hegv.fypp"
#:include "src/mfi/lapack/heevd.fypp"
#:include "src/mfi/lapack/potrf_potri.fypp"
#:include "src/mfi/lapack/potrs.fypp"
#:include "src/mfi/lapack/pocon.fypp"
#:set COLLECT = [                            &
    ('?geqrf',  DEFAULT_TYPES, geqrf_gerqf), &
    ('?gerqf',  DEFAULT_TYPES, geqrf_gerqf), &
    ('?getrf',  DEFAULT_TYPES, getrf),       &
    ('?getri',  DEFAULT_TYPES, getri),       &
    ('?getrs',  DEFAULT_TYPES, getrs),       &
    ('?hetrf',  COMPLEX_TYPES, hetrf),       &
    ('?hegv',   COMPLEX_TYPES, hegv),        &
    ('?heevd',  COMPLEX_TYPES, heevd),       &
    ('?gesvd',  DEFAULT_TYPES, gesvd),       &
    ('?potrf',  DEFAULT_TYPES, potrf_potri), &
    ('?potri',  DEFAULT_TYPES, potrf_potri), &
    ('?potrs',  DEFAULT_TYPES, potrs),       &
    ('?pocon',  DEFAULT_TYPES, pocon),       &
]
#:endmute
!> Modern fortran interfaces for LAPACK
module mfi_lapack
use iso_fortran_env
use f77_lapack
use f77_lapack, only: mfi_lartg => f77_lartg
implicit none

#:for name, supported_types, code in COLLECT
$:mfi_interface(name, supported_types)
#:endfor

contains

#:for name, supported_types, code in COLLECT
$:mfi_implement(name, supported_types, code)
#:endfor

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module



================================================
FILE: src/mfi/blas/asum_nrm2.fypp
================================================
#:def asum_nrm2(MFI_NAME,F77_NAME,pfxs)
#:set A, B = get_types(pfxs)
pure function ${MFI_NAME}$(x, incx)
    ${type(A)}$ :: ${MFI_NAME}$
@:args(${type(B)}$, in, x(:))
@:optional(integer, in, incx)
    integer :: n
@:defaults(incx=1)
    n = size(x)
    ${MFI_NAME}$ = ${F77_NAME}$(n, x, local_incx)
end function
#:enddef



================================================
FILE: src/mfi/blas/axpy.fypp
================================================
#:def axpy(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(x, y, a, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:))
@:args(${type(wp)}$, inout, y(:))
@:optional(${type(wp)}$, in, a)
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(a=1.0_wp, incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,local_a,x,local_incx,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/copy_swap.fypp
================================================
#:def copy_swap(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(x, y, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:))
@:args(${type(wp)}$, inout, y(:))
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,x,local_incx,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/dot_product.fypp
================================================
#:def dot_product(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure function ${MFI_NAME}$(x, y, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
    ${type(wp)}$ :: ${MFI_NAME}$
@:args(${type(wp)}$, in, x(:), y(:))
    integer :: n
@:optional(integer, in, incx, incy)
@:defaults(incx=1, incy=1)
    N = size(X)
    ${MFI_NAME}$ = ${F77_NAME}$(n,x,local_incx,y,local_incy)
end function
#:enddef



================================================
FILE: src/mfi/blas/gbmv.fypp
================================================
#:def gbmv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, y, kl, m, alpha, beta, trans, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:), x(:))
@:args(${type(wp)}$, inout, y(:))
@:optional(character, in, trans)
@:optional(${type(wp)}$,  in, alpha, beta)
@:optional(integer,   in, kl, m, incx,  incy)
    integer :: n, ku, lda
    n = size(a,2)
    lda = max(1,size(a,1))
@:defaults(kl=(lda-1)/2, m=n, trans='N', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    ku = lda-local_kl-1
    call ${F77_NAME}$(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/gemm.fypp
================================================
#:def gemm(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, c, transa, transb, alpha, beta)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:), b(:,:))
@:args(${type(wp)}$, inout, c(:,:))
@:optional(character, in, transa, transb)
@:optional(${type(wp)}$,  in, alpha, beta)
    integer :: m, n, k, lda, ldb, ldc
@:defaults(transa='N', transb='N', alpha=1.0_wp, beta=0.0_wp)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call ${F77_NAME}$(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/gemv.fypp
================================================
#:def gemv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, y, trans, alpha, beta, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:), x(:))
@:args(${type(wp)}$, inout, y(:))
@:optional(character, in, trans)
@:optional(${type(wp)}$,  in, alpha, beta)
@:optional(integer,   in, incx,  incy)
    integer :: m, n, lda
@:defaults(trans='N', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call ${F77_NAME}$(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/ger_gerc_geru.fypp
================================================
#:def ger_gerc_geru(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, y, alpha, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:), y(:))
@:args(${type(wp)}$, inout, a(:,:))
@:optional(${type(wp)}$,  in, alpha)
@:optional(integer,   in, incx, incy)
    integer :: m, n, lda
@:defaults(alpha=1.0_wp, incx=1, incy=1)
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call ${F77_NAME}$(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/hbmv_sbmv.fypp
================================================
#:def hbmv_sbmv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, y, uplo, alpha, beta, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in, x(:), a(:,:))
@:args(${type(wp)}$, inout, y(:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha, beta)
@:optional(integer,   in, incx, incy)
    integer :: n, k, lda
@:defaults(uplo='U', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/hemm_symm.fypp
================================================
#:def hemm_symm(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, c, side, uplo, alpha, beta)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:), b(:,:))
@:args(${type(wp)}$, inout, c(:,:))
@:optional(character, in, side,  uplo)
@:optional(${type(wp)}$,  in, alpha, beta)
    integer :: m, n, lda, ldb, ldc
@:defaults(side='L', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call ${F77_NAME}$(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/hemv_symv.fypp
================================================
#:def hemv_symv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, y, uplo, alpha, beta, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in, x(:), a(:,:))
@:args(${type(wp)}$, inout, y(:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha, beta)
@:optional(integer,   in, incx, incy)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/her.fypp
================================================
#:def her(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, uplo, alpha, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:))
@:args(${type(wp)}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(real(wp),  in, alpha)
@:optional(integer,   in, incx)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/her2k.fypp
================================================
#:def her2k(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:))
@:args(${type(wp)}$, in,    b(:,:))
@:args(${type(wp)}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(${type(wp)}$,  in, alpha)
@:optional(real(wp),  in, beta)
    integer :: n, k, lda, ldb, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/her_syr2.fypp
================================================
#:def her_syr2(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, y, uplo, alpha, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:), y(:))
@:args(${type(wp)}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha)
@:optional(integer,   in, incx, incy)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, incx=1, incy=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/herk.fypp
================================================
#:def herk(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:))
@:args(${type(wp)}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(real(wp),  in, alpha, beta)
    integer :: n, k, lda, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/hpmv_spmv.fypp
================================================
#:def hpmv_spmv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(ap, x, y, uplo, alpha, beta, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:), ap(:))
@:args(${type(wp)}$, inout, y(:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha, beta)
@:optional(integer,   in, incx, incy)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/hpr.fypp
================================================
#:def hpr(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(ap, x, uplo, alpha, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:))
@:args(${type(wp)}$, inout, ap(:))
@:optional(character, in, uplo)
@:optional(real(wp),  in, alpha)
@:optional(integer,   in, incx)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/hpr_spr2.fypp
================================================
#:def hpr_spr2(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(ap, x, y, uplo, alpha, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:), y(:))
@:args(${type(wp)}$, inout, ap(:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha)
@:optional(integer,   in, incx, incy)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, incx=1, incy=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/iamin_iamax.fypp
================================================
#:def iamin_iamax(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure function ${MFI_NAME}$(x, incx)
@:parameter(integer, wp=${kind(wp)}$)
    integer :: ${MFI_NAME}$
@:args(${type(wp)}$, in, x(:))
@:optional(integer, in, incx)
    integer :: n
@:defaults(incx=1)
    n = size(x)
    ${MFI_NAME}$ = ${F77_NAME}$(n,x,local_incx)
end function
#:enddef



================================================
FILE: src/mfi/blas/lamch.fypp
================================================
#:def lamch(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure function ${MFI_NAME}$(cmach, kind) result(res)
@:parameter(integer, wp=${kind(wp)}$)
@:args(character, in, cmach)
@:args(${type(wp)}$, in, kind)
    !! Just a kind placeholder
    ${type(wp)}$ :: res
    res = ${F77_NAME}$(cmach)
end function
#:enddef



================================================
FILE: src/mfi/blas/rot.fypp
================================================
#:def rot(MFI_NAME,F77_NAME,pfxs)
#:set A, B = get_types(pfxs)
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
#:if type(A) == real(A)
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
#:elif type(A) == complex(A)
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
#:endif
!>```
pure subroutine ${MFI_NAME}$(x, y, c, s, incx, incy)
@:parameter(integer, wp=${kind(A)}$)
@:args(${type(A)}$, inout, x(:), y(:))
@:args(${real(A)}$, in, c)
@:args(${type(B)}$, in, s)
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(incx=1, incy=1)
    n = size(x)
    call ${F77_NAME}$(n,x,local_incx,y,local_incy,c,s)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/rotm.fypp
================================================
#:def rotm(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(x, y, param, incx, incy)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, x(:), y(:))
@:args(${type(wp)}$, in, param(5))
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,x,local_incx,y,local_incy,param)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/scal.fypp
================================================
#:def scal(MFI_NAME,F77_NAME,pfxs)
#:set A, B = get_types(pfxs)
!> ${MFI_NAME.upper()}$ scales a vector by a constant.
pure subroutine ${MFI_NAME}$(a, x, incx)
@:args(${type(A)}$, inout, x(:))
@:args(${type(B)}$, in,    a)
@:optional(integer, in, incx)
    integer :: n
@:defaults(incx=1)
    n = size(x)
    call ${F77_NAME}$(n,a,x,local_incx)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/spr.fypp
================================================
#:def spr(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(ap, x, uplo, alpha, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:))
@:args(${type(wp)}$, inout, ap(:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha)
@:optional(integer,   in, incx)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/syr.fypp
================================================
#:def syr(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, uplo, alpha, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    x(:))
@:args(${type(wp)}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(${type(wp)}$,  in, alpha)
@:optional(integer,   in, incx)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/syr2k.fypp
================================================
#:def syr2k(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:))
@:args(${type(wp)}$, in,    b(:,:))
@:args(${type(wp)}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(${type(wp)}$,  in, alpha, beta)
    integer :: n, k, lda, ldb, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/syrk.fypp
================================================
#:def syrk(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:))
@:args(${type(wp)}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(${type(wp)}$,  in, alpha, beta)
    integer :: n, k, lda, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/tbmv_tbsv.fypp
================================================
#:def tbmv_tbsv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, uplo, trans, diag, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in, a(:,:))
@:args(${type(wp)}$, inout, x(:))
@:optional(character, in, uplo, trans, diag)
@:optional(integer,   in, incx)
    integer :: n, k, lda
@:defaults(uplo='U', trans='N', diag='N', incx=1)
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/tpmv_tpsv.fypp
================================================
#:def tpmv_tpsv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(ap, x, uplo, trans, diag, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    ap(:))
@:args(${type(wp)}$, inout, x(:))
@:optional(character, in, uplo, trans, diag)
@:optional(integer,   in, incx)
    integer :: n
@:defaults(uplo='U', trans='N', diag='N', incx=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/trmm_trsm.fypp
================================================
#:def trmm_trsm(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, side, uplo, transa, diag, alpha)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in,    a(:,:))
@:args(${type(wp)}$, inout, b(:,:))
@:optional(character, in, side, uplo, transa, diag)
@:optional(${type(wp)}$,  in, alpha)
    integer :: m, n, lda, ldb
@:defaults(side='L', uplo='U', transa='N', diag='N', alpha=1.0_wp)
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call ${F77_NAME}$(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
#:enddef



================================================
FILE: src/mfi/blas/trmv_trsv.fypp
================================================
#:def trmv_trsv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, x, uplo, trans, diag, incx)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, in, a(:,:))
@:args(${type(wp)}$, inout, x(:))
@:optional(character, in, uplo, trans, diag)
@:optional(integer,   in, incx)
    integer :: n, lda
@:defaults(uplo='U', trans='N', diag='N', incx=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/geqrf_gerqf.fypp
================================================
#:def geqrf_gerqf(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, tau, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(:,:))
    ${type(wp)}$, intent(out), optional, target :: tau(:)
@:optional(integer, out, info)
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    ${type(wp)}$, pointer :: local_tau(:), work(:)
    ${type(wp)}$, target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call ${F77_NAME}$(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call ${F77_NAME}$(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/gesvd.fypp
================================================
#:def gesvd(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, s, u, vt, ww, job, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(:,:))
@:args(${real(wp)}$, out, s(:))
    ${type(wp)}$,      intent(out), optional, target :: u(:,:), vt(:,:)
    ${real(wp)}$, intent(out), optional, target :: ww(:)
@:optional(character, in, job)
@:optional(integer, out, info)
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    ${type(wp)}$,      target  :: s_work(1), l_a2(1,1)
    ${type(wp)}$,      pointer :: local_u(:,:), local_vt(:,:), work(:)
#:if type(wp) == complex(wp)
    ${real(wp)}$, pointer :: rwork(:)
#:endif
@:defaults(job='N')
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
@:optval(present(u),  ldu,  max(1,size(u,1)),  1)
@:optval(present(vt), ldvt, max(1,size(vt,1)), 1)
    if (present(u)) then
        if (size(u,2) == m) then
            jobu = 'A'
        else
            jobu = 'S'
        end if
        local_u => u
    else
        if (local_job == 'u' .or. local_job == 'U') then
            jobu = 'O'
        else
            jobu = 'N'
        end if
        local_u => l_a2
    end if
    if (present(vt)) then
        if (size(vt,1) == n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        end if
        local_vt => vt
    else
        if (local_job == 'v' .or. local_job == 'V') then
            jobvt = 'O'
        else
            jobvt = 'N'
        end if
        local_vt => l_a2
    end if
    allocation_status = 0
    lwork = -1
#:if type(wp) == complex(wp)
    allocate(rwork(5*min(m,n)), stat=allocation_status)
    call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,rwork,local_info)
#:else
    call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,local_info)
#:endif
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
#:if type(wp) == complex(wp)
        call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,rwork,local_info)
#:else
        call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,local_info)
#:endif
    else
        local_info = -1000
    end if

    if (present(ww)) then
        ww = real(work(2:min(m,n)-1))
    end if
    deallocate(work, stat=deallocation_status)
404 continue
#:if type(wp) == complex(wp)
    deallocate(rwork, stat=deallocation_status)
#:endif
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/getrf.fypp
================================================
#:def getrf(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, ipiv, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(:,:))
    integer, intent(out), optional, target :: ipiv(:)
@:optional(integer, out, info)
    integer :: m, n, lda, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(min(m,n)), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call ${F77_NAME}$(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/getri.fypp
================================================
#:def getri(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, ipiv, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(:,:))
@:args(integer,     in, ipiv(:))
    ${type(wp)}$, pointer :: work(:)
    ${type(wp)}$ :: s_work(1)
@:optional(integer, out, info)
    integer :: n, lda, lwork, allocation_status, deallocation_status
    lda = max(1,size(a,1))
    n = size(a,2)
    lwork = -1
    call ${F77_NAME}$(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call ${F77_NAME}$(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$',-local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/getrs.fypp
================================================
#:def getrs(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a,ipiv,b,trans,info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(:,:))
@:args(${type(wp)}$, inout, b(:,:))
@:args(integer,     in, ipiv(:))
@:optional(integer, out, info)
@:optional(character, in, trans)
    integer :: n, nrhs, lda, ldb
@:defaults(trans='N')
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call ${F77_NAME}$(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$',-local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/heevd.fypp
================================================
#:def heevd(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, w, jobz, uplo, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(:,:))
@:args(${real(wp)}$, out, w(:))
@:optional(integer, out, info)
@:optional(character, in, jobz, uplo)
    ${type(wp)}$,      pointer :: work(:)
    ${real(wp)}$, pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    ${type(wp)}$      :: s_work(1)
    ${real(wp)}$ :: s_rwork(1)
    integer       :: s_iwork(1)
    integer :: n, lda, lwork, lrwork, liwork, allocation_status, deallocation_status
@:defaults(jobz='N', uplo='U')
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1

    call ${F77_NAME}$(local_jobz,local_uplo,n,a,lda,w, &
                      s_work,lwork,s_rwork,lrwork,s_iwork,liwork,local_info)
    if (local_info /= 0) goto 404
    lwork  = int(s_work(1))
    lrwork = int(s_rwork(1))
    liwork = int(s_iwork(1))

    allocate(iwork(liwork), stat=allocation_status)

    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        allocate(work(lwork),   stat=allocation_status)
        call ${F77_NAME}$(local_jobz,local_uplo,n,a,lda,w, &
                      work,lwork,rwork,lrwork,iwork,liwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(iwork, stat=deallocation_status)
    deallocate(rwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/heevr.fypp
================================================
#:mute
#:def heevr(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, w, uplo, z, vl, vu, il, iu, m, isuppz, abstol, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,      inout, a(:,:))
@:args(${real(wp)}$, out,   w(:))
@:optional(character,     in,    uplo)
@:optional(${real(wp)}$, in,    vl, vu, abstol)
@:optional(integer,       in,    il, iu)
@:optional(integer,       out,   m)
@:optional(integer,       out,   info)
    ${type(wp)}$,      intent(out), optional, target :: z(:,:)
    integer,       intent(out), optional, target :: isuppz(:)
    integer,       pointer :: local_isuppz(:)
    ${type(wp)}$,      pointer :: local_z(:,:)
    ${type(wp)}$,      pointer :: work(:)
    ${real(wp)}$, pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    ${type(wp)}$      :: s_work(1)
    ${real(wp)}$ :: s_rwork(1)
    integer       :: s_iwork(1)
    character     :: jobz, range
    integer :: n, lda, ldz, lwork, lrwork, liwork, allocation_status, deallocation_status
    integer, target  :: dummy_rank_1(1)
    ${type(wp)}$, target :: dummy_rank_2(1,1)

    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1

@:defaults(uplo='U', abstol=0.0_wp, il=1, iu=N, vl=huge(vl), vu=huge(vu))

    if (present(Z)) then
        jobz    = 'V'
        ldz     = max(1,size(z,1))
        local_z = Z
        if (present(isuppz)) then
            local_isuppz => isuppz
        else
            allocate(local_isuppz(2*N), stat=allocation_status)
        end if
    else
        jobz    = 'N'
        ldz     = 1
        local_z = dummy_rank_2
        if (present(isuppz)) then
            local_info = -1001
            goto 404
        else
            local_isuppz => dummy_rank_1
        end if
    end if

    if (present(vl) .or. present(vu) &
  .and. present(il) .or. present(iu)) then
        local_info = -1001
        goto 404
    else if (present(vl) .or. present(vu)) then
        range = 'V'
    else if (present(il) .or. present(iu)) then
        range = 'I'
    else
        range = 'A'
    end if

    call ${F77_NAME}$(jobz,range,local_uplo,n,a,lda, &
                      local_vl,local_vu, &
                      local_il,local_iu, &
                      local_abstol,local_m,w,local_z,ldz,local_isuppz, &
                      s_work,lwork,s_rwork,lrwork, &
                      s_iwork,liwork,local_info)
    if (local_info /= 0) goto 405
    lwork  = int(s_work(1))
    lrwork = int(s_rwork(1))
    liwork = int(s_iwork(1))

    if (allocation_status == 0) allocate(iwork(liwork), stat=allocation_status)
    if (allocation_status == 0) allocate(rwork(lrwork), stat=allocation_status)
    if (allocation_status == 0) allocate(work(lwork),   stat=allocation_status)
    if (allocation_status == 0) then
        call ${F77_NAME}$(jobz,range,local_uplo,n,a,lda, &
                          local_vl,local_vu, &
                          local_il,local_iu, &
                          local_abstol,local_m,w,local_z,ldz,local_isuppz, &
                          work,lwork,rwork,lrwork, &
                          iwork,liwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(iwork, stat=deallocation_status)
    deallocate(rwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)
405 continue
    if (present(z) .and. .not. present(isuppz)) then
        deallocate(local_isuppz, stat=deallocation_status)
    end if
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef
#:endmute



================================================
FILE: src/mfi/lapack/hegv.fypp
================================================
#:def hegv(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, w, itype, jobz, uplo, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(:,:), b(:,:))
@:args(${real(wp)}$, out, w(:))
@:optional(integer,   in, itype)
@:optional(character, in, jobz, uplo)
@:optional(integer, out, info)
    ${type(wp)}$,      pointer :: work(:)
    ${real(wp)}$, pointer :: rwork(:)
    ${type(wp)}$ :: s_work(1)
    integer :: n, lda, ldb, lwork, allocation_status, deallocation_status
@:defaults(itype=1, jobz='N', uplo='U')
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    allocation_status = 0
    allocate(rwork(max(1,3*N-2)), stat=allocation_status)
    lwork = -1
    call ${F77_NAME}$(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call ${F77_NAME}$(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/hetrf.fypp
================================================
#:def hetrf(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, uplo, ipiv, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,  inout, a(:,:))
    integer, intent(out), optional, target :: ipiv(:)
    integer, pointer :: local_ipiv(:)
@:optional(character, in, uplo)
@:optional(integer, out, info)
    integer   :: n, lda, lwork, allocation_status, deallocation_status
    ${type(wp)}$, target :: s_work(1)
    ${type(wp)}$, pointer :: work(:)
@:defaults(uplo='U')
    lda = max(1,size(a,1))
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(n), stat=allocation_status)
    end if
    lwork = -1
    call ${F77_NAME}$(local_uplo,n,a,lda,local_ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (.not. present(ipiv)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$',-local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/pocon.fypp
================================================
#:mute

!subroutine mfi_pocon (
!   type(wp), dimension(:,:)     a,
!   real(wp)                 anorm,
!   real(wp)                 rcond,
!   character                 uplo,
!   integer                   info
!)

#:def pocon(MFI_NAME, F77_NAME, pfxs)
#:set wp = pfxs[0]
!> Estimates the reciprocal of the condition number of a real symmetric / complex Hermitian positive definite matrix using the Cholesky factorization computed by ?POTRF
pure subroutine ${MFI_NAME}$(a, anorm, rcond, uplo, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(:,:))
@:args(${real(wp)}$, in, anorm)
@:args(${real(wp)}$, out, rcond)
@:optional(character, in, uplo)
@:optional(integer, out, info)
    integer :: n, lda, allocation_status, deallocation_status
    ${type(wp)}$, pointer :: work(:)
#:if type(wp) == real(wp)
    integer, pointer :: xwork(:)
#:elif type(wp) == complex(wp)
    ${real(wp)}$, pointer :: xwork(:)
#:endif
@:defaults(uplo='U')
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    allocate(xwork(n), stat=allocation_status)
    if (allocation_status == 0) allocate(work(3*n), stat=allocation_status)

    if (allocation_status == 0) then
        call ${F77_NAME}$(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$',-local_info)
    end if
end subroutine
#:enddef
#:endmute



================================================
FILE: src/mfi/lapack/potrf_potri.fypp
================================================
#:def potrf_potri(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, info, uplo)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(integer, out, info)
    integer :: n, lda
@:defaults(uplo='U')
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('${F77_NAME}$', local_info)
    end if
end subroutine
#:enddef



================================================
FILE: src/mfi/lapack/potrs.fypp
================================================
#:def potrs(MFI_NAME,F77_NAME,pfxs)
#:set wp = pfxs[0]
pure subroutine ${MFI_NAME}$(a, b, uplo, info)
@:parameter(integer, wp=${kind(wp)}$)
@:args(${type(wp)}$,    in, a(:,:))
@:args(${type(wp)}$, inout, b(:,:))
@:optional(character, in, uplo)
@:optional(integer, out, info)
    integer :: n, nrhs, lda, ldb
@:defaults(uplo='U')
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call ${F77_NAME}$(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$',-local_info)
    end if
end subroutine
#:enddef



================================================
FILE: test/blas.fpp
================================================
#:mute
#:include "common.fpp"
#:include "test/blas/asum_nrm2.fypp"
#:include "test/blas/iamin_iamax.fypp"
#:include "test/blas/dot_product.fypp"
#:include "test/blas/axpy.fypp"
#:include "test/blas/rot.fypp"
#:include "test/blas/rotg.fypp"
#:include "test/blas/copy.fypp"
#:include "test/blas/swap.fypp"
#:include "test/blas/scal.fypp"
#:include "test/blas/lamch.fypp"
#:include "test/blas/gemv.fypp"
#:include "test/blas/gemm.fypp"
#:set COLLECT = [                                             &
    ('?lamch',REAL_TYPES,                       lamch),       &
    ('?dot',  REAL_TYPES,                       dot_product), &
    ('?dotc', COMPLEX_TYPES,                    dot_product), &
    ('?dotu', COMPLEX_TYPES,                    dot_product), &
    ('?copy', DEFAULT_TYPES,                    copy),        &
    ('?swap', DEFAULT_TYPES,                    swap),        &
    ('?axpy', DEFAULT_TYPES,                    axpy),        &
    ('?gemv', DEFAULT_TYPES,                    gemv),        &
    ('?gemm', DEFAULT_TYPES,                    gemm),        &
    ('?asum', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),   &
    ('?nrm2', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),   &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),         &
    ('?rotg', DEFAULT_TYPES, rotg),                           &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),        &
]
#:endmute
program main
use iso_fortran_env
implicit none
#:for name, supported_types, code in COLLECT
$:test_run(name, supported_types)
#:endfor

#:if defined('MFI_EXTENSIONS')
$:test_run('i?amin',DEFAULT_TYPES)
$:test_run('i?amax',DEFAULT_TYPES)
#:endif

contains

#:for name, supported_types, code in COLLECT
$:test_implement(name, supported_types,code)
#:endfor

#:if defined('MFI_EXTENSIONS')
$:test_implement('i?amin',DEFAULT_TYPES, iamin_iamax)
$:test_implement('i?amax',DEFAULT_TYPES, iamin_iamax)
#:endif

pure subroutine assert(test, msg)
    logical, intent(in) :: test
    character(*), intent(in) :: msg
    if (.not. test) then
        error stop msg
    end if
end subroutine

end program



================================================
FILE: test/lapack.fpp
================================================
#:include "common.fpp"

program test_mfi_lapack
    use iso_fortran_env
    use mfi_lapack
    use f77_lapack
    implicit none
    integer, parameter :: N = 2000
    real(REAL64) :: S(N,N)
    integer :: i, j, info

    call test_geqrf
    call test_gesvd
    call test_potrf
    call test_potri

contains

    subroutine symmetric
        call random_number(S)
        do j=1,N
            do i=j+1,N
                S(i,j) = S(j,i)
            end do
        end do
    end subroutine

    subroutine positive_definite
        call symmetric
        S = matmul(S,S)
    end subroutine

    subroutine test_gesvd
        real(REAL64) :: A(3,5)
        real(REAL64) :: S(3), ES(3)
        ES = [real(REAL64) :: 14.0828, 10.1124, 3.92600]
        A(1,:) = [real(REAL64) :: -2,  1, -5, 11, 2]
        A(2,:) = [real(REAL64) ::  2, -2,  0,  1, 4]
        A(3,:) = [real(REAL64) ::  6, -8,  0,  6, 0]
        @:timeit("mfi_gesvd: ", { call mfi_gesvd(A,S) })
        call assert(all(abs(S-ES) < 1e-4))
    end subroutine

    subroutine test_geqrf
        integer, parameter :: wp=REAL64
        integer, parameter :: N=6, M=2
        real(wp) :: A(N,M), B(N,M)
        real(wp) :: tau(min(N,M)), tau_(min(N,M))

        A(1,:) = [  .000000_wp,  2.000000_wp]
        A(2,:) = [ 2.000000_wp, -1.000000_wp]
        A(3,:) = [ 2.000000_wp, -1.000000_wp]
        A(4,:) = [  .000000_wp,  1.500000_wp]
        A(5,:) = [ 2.000000_wp, -1.000000_wp]
        A(6,:) = [ 2.000000_wp, -1.000000_wp]

        B(1,:) = [ -4.000000_wp, 2.000000_wp]
        B(2,:) = [   .500000_wp, 2.500000_wp]
        B(3,:) = [   .500000_wp,  .285714_wp]
        B(4,:) = [   .000000_wp, -.428571_wp]
        B(5,:) = [   .500000_wp,  .285714_wp]
        B(6,:) = [   .500000_wp,  .285714_wp]

        tau_ = [1.0_wp, 1.4_wp]

        @:timeit("mfi_geqrf: ", { call mfi_geqrf(A, tau) })
        call assert(all(abs(A-B) < 1e-6))
        call assert(all(abs(tau-tau_) < 1e-6))
    end subroutine

    subroutine test_potrf
        call positive_definite
        @:timeit("f77_potrf: ", { call f77_potrf('U',N,S,N,info) })
        call assert(info == 0)

        call positive_definite
        @:timeit("mfi_potrf: ", { call mfi_potrf(S,info)         })
        call assert(info == 0)
    end subroutine

    subroutine test_potri
        call positive_definite
        call f77_potrf('U',N,S,N,info)
        @:timeit("f77_potri: ", { call f77_potri('U',N,S,N,info) })
        call assert(info == 0)

        call positive_definite
        call mfi_potrf(S,info)
        @:timeit("mfi_potri: ", { call mfi_potri(S,info)         })
        call assert(info == 0)
    end subroutine

    pure subroutine assert(test)
        logical, intent(in) :: test
        if (.not. test) then
            error stop 'assertion failed'
        end if
    end subroutine

end program



================================================
FILE: test/blas/asum_nrm2.fypp
================================================
#:def asum_nrm2(f77,f90,mfi,pfxs)
#:set A, B = get_types(pfxs)
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(A)}$
    integer, parameter :: N = 20
    ${type(B)}$ :: array(N)
    ${type(A)}$ :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = ${f77}$(N, array, 1)
    res(3) = ${f90}$(N, array, 1)
    res(2) = mfi_${f77}$(array)
    res(4) = ${mfi}$(array)
    call assert(all(res == res(1)), "different results for sequential array")

    $:random_number(type(B),'array','(N)')
    res(1) = ${f77}$(N, array, 1)
    res(2) = ${f90}$(N, array, 1)
    res(3) = mfi_${f77}$(array)
    res(4) = ${mfi}$(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
#:enddef



================================================
FILE: test/blas/axpy.fypp
================================================
#:def axpy(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20
    ${real(wp)}$ :: rnd_vector(N), rnd
    ${type(wp)}$ :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    ${type(wp)}$ :: alpha

#:if type(wp) == complex(wp)
    call random_number(rnd_vector)
    X%re = rnd_vector
    call random_number(rnd_vector)
    X%im = rnd_vector
    call random_number(rnd_vector)
    Y%re = rnd_vector
    call random_number(rnd_vector)
    Y%im = rnd_vector
    call random_number(rnd)
    alpha%re = rnd
    call random_number(rnd)
    alpha%im= rnd
#:else
    call random_number(X)
    call random_number(Y)
    call random_number(alpha)
#:endif

    x_in = X
    y_in = Y
    call ${f77}$(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call ${f90}$(N, alpha, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_${f77}$(x_in,y_in,alpha)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call ${mfi}$(x_in, y_in, alpha)

    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
#:enddef



================================================
FILE: test/blas/copy.fypp
================================================
#:mute
#:def copy(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20

    ${type(wp)}$ :: x(N),    y(N),    &
                    x_in(N), y_in(N), &
                    x_rf(N), y_rf(N)

    $:random_number(type(wp),'x','(N)')
    $:random_number(type(wp),'y','(N)')

    x_in = x
    y_in = y
    ! The test is always against the original
    call ${f77}$(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call ${f90}$(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_${f77}$(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call ${mfi}$(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
#:enddef
#:endmute



================================================
FILE: test/blas/dot_product.fypp
================================================
#:mute
#:def dot_product(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20

    ${type(wp)}$ :: res, ref

    ${type(wp)}$ :: x(N), y(N)

    $:random_number(type(wp),'X','(N)')
    $:random_number(type(wp),'Y','(N)')

    ! The test is always against the original
    ref = ${f77}$(N, x, 1, y, 1)

    res = ${f90}$(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_${f77}$(x, y)
    call assert(ref == res, "different results")

    res = ${mfi}$(x, y)
    call assert(ref == res, "different results")

end subroutine
#:enddef
#:endmute



================================================
FILE: test/blas/gemm.fypp
================================================
#:def gemm(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    ${type(wp)}$ :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    ${type(wp)}$ :: alpha, beta
    character :: transa, transb
    integer :: i, j

    $:random_number(type(wp),'A','(N,N)')
    $:random_number(type(wp),'B','(N,N)')
    $:random_number(type(wp),'C','(N,N)')
    $:random_number(type(wp),'alpha')
    $:random_number(type(wp),'beta')

    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call ${f77}$(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call ${f90}$(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_${f77}$(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call ${mfi}$(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
#:enddef



================================================
FILE: test/blas/gemv.fypp
================================================
#:def gemv(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20
    ${type(wp)}$ :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    ${type(wp)}$ :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

    $:random_number(type(wp),'M','(N,N)')
    $:random_number(type(wp),'X','(N)')
    $:random_number(type(wp),'Y','(N)')
    $:random_number(type(wp),'alpha')
    $:random_number(type(wp),'beta')


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call ${f77}$(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call ${f90}$(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_${f77}$(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call ${mfi}$(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
#:enddef



================================================
FILE: test/blas/iamin_iamax.fypp
================================================
#:def iamin_iamax(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20
    ${real(wp)}$ :: rnd(N)
    ${type(wp)}$ :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = ${f77}$(N, array, 1)
    res(2) = mfi_${f77}$(array)
    res(3) = ${f90}$(N, array, 1)
    res(4) = ${mfi}$(array)
    call assert(all(res == res(1)), "different results for sequential array")

#:if type(wp) == complex(wp)
    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
#:else
    call random_number(array)
#:endif
    res(1) = ${f77}$(N, array, 1)
    res(2) = ${f90}$(N, array, 1)
    res(3) = mfi_${f77}$(array)
    res(4) = ${mfi}$(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
#:enddef



================================================
FILE: test/blas/lamch.fypp
================================================
#:mute
#:def lamch(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$
    use mfi_blas, only: ${mfi}$

    integer, parameter :: wp = ${kind(wp)}$

    integer, parameter :: N = 20
    character, parameter :: options(*) = ['E','e', &
                                          'S','s', &
                                          'B','b', &
                                          'P','p', &
                                          'N','n', &
                                          'R','r', &
                                          'M','m', &
                                          'U','u', &
                                          'L','l', &
                                          'O','o']
    ${type(wp)}$ :: a, b
    integer :: i

    do i=1,size(options)
        a = ${f77}$(options(i))
        b = ${mfi}$(options(i),1.0_wp)
        call assert(a == b, "different results for option "//options(i))
    end do

end subroutine
#:enddef
#:endmute



================================================
FILE: test/blas/rot.fypp
================================================
#:mute
#:def rot(f77,f90,mfi,pfxs)
#:set A, B = get_types(pfxs)
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer,     parameter :: wp = ${kind(A)}$
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    ${type(A)}$ :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    ${real(A)}$ :: c
    ${type(B)}$ :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    $:random_number(type(A),'X','(N)')
    $:random_number(type(A),'Y','(N)')

    c = cos(angle)
    s = sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call ${f77}$(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call ${f90}$(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_${f77}$(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call ${mfi}$(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
#:enddef
#:endmute



================================================
FILE: test/blas/rotg.fypp
================================================
#:def rotg(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$
    use mfi_blas, only: ${mfi}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 200
    ${type(wp)}$ :: a, b, s
    ${real(wp)}$ :: c

    ${type(wp)}$ :: a_in, b_in, s_in
    ${real(wp)}$ :: c_in

    ${type(wp)}$ :: a_rf, b_rf, s_rf
    ${real(wp)}$ :: c_rf
    integer :: i

    $:random_number(type(wp),'a')
    $:random_number(type(wp),'b')
    $:random_number(real(wp),'c')
    $:random_number(type(wp),'s')

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call ${f77}$(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call ${mfi}$(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
#:enddef



================================================
FILE: test/blas/scal.fypp
================================================
#:mute
#:def scal(f77,f90,mfi,pfxs)
#:set A, B = get_types(pfxs)
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer,     parameter :: wp = ${kind(A)}$
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    ${type(A)}$ :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    ${type(B)}$ :: alpha

    $:random_number(type(A),'X','(N)')
    $:random_number(type(B),'alpha')

    ! The test is always against the original
    x_in = x
    call ${f77}$(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call ${f90}$(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_${f77}$(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call ${mfi}$(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine
#:enddef
#:endmute



================================================
FILE: test/blas/swap.fypp
================================================
#:mute
#:def swap(f77,f90,mfi,pfxs)
#:set wp = pfxs[0]
subroutine test_${f77}$
    use f77_blas, only: ${f77}$, ${f90}$
    use mfi_blas, only: ${mfi}$, mfi_${f77}$

    integer, parameter :: wp = ${kind(wp)}$
    integer, parameter :: N = 20

    ${real(wp)}$ :: rnd(N)

    ${type(wp)}$ :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

#:if type(wp) == complex(wp)
    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd
#:else
    call random_number(X)
    call random_number(Y)
#:endif

    x_in = x
    y_in = y
    ! The test is always against the original
    call ${f77}$(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call ${f90}$(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_${f77}$(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call ${mfi}$(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
#:enddef
#:endmute



================================================
FILE: .github/workflows/docs-deployment.yml
================================================
name: docs-deployment

on: [push]

jobs:
  Build:
    runs-on: ubuntu-latest

    steps:
    - name: Checkout 🛎️
      uses: actions/checkout@v2.3.1

    - name: Set up Python 3.x
      uses: actions/setup-python@v1
      with:
        python-version: 3.x

    - name: Install fypp
      run: pip install --upgrade fypp

    - name: Install ford
      run: pip install --upgrade ford

    - name: Generate fpm package 🔧
      run: |
        make FYPPFLAGS=-DMFI_EXTENSIONS
        ford ford.md

    - name: Deploy docs 🚀
      uses: JamesIves/github-pages-deploy-action@4.1.5
      if: github.event_name != 'pull_request'
      with:
        BRANCH: docs
        FOLDER: api-reference



================================================
FILE: .github/workflows/fpm-deployment.yml
================================================
name: fpm-deployment

on: [push, pull_request]
env:
  PROJECT_NAME: "mfi"

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        include:
          - os: ubuntu-latest
            toolchain: {compiler: gcc, version: 13}
    steps:
    - name: Checkout 🛎️
      uses: actions/checkout@v2.3.1

    - name: Set up Python 3.x
      uses: actions/setup-python@v1
      with:
        python-version: 3.x

    - name: Install fypp
      run: pip install --upgrade fypp

    - name: Generate ${{ env.PROJECT_NAME }}-fpm package 🔧
      run: |
        make FYPPFLAGS=-DMFI_EXTENSIONS
        mkdir mfi-fpm
        cp -R src   mfi-fpm
        cp -R test  mfi-fpm
        cp fpm.toml mfi-fpm
        cp LICENSE  mfi-fpm
        find mfi-fpm/src  -type f ! -name "*.f90" -delete
        find mfi-fpm/test -type f ! -name "*.f90" -delete
    
    - name: setup fortran
      uses: fortran-lang/setup-fortran@main
      id: setup-fortran
      with:
        compiler: ${{ matrix.toolchain.compiler }}
        version: ${{ matrix.toolchain.version }}

    - name: Install Lapack and Blas
      run: |
        sudo apt install libblas-dev liblapack-dev libtmglib-dev

    - name: Install fpm latest release
      uses: fortran-lang/setup-fpm@v5
      with:
        fpm-version: 'v0.10.0'

    - name: Run fpm build ⚙
      run: |
        fpm build --profile release

    # Update and deploy the f90 files generated by github-ci to the `PROJECT_NAME-fpm` branch.
    - name: Deploy 🚀
      uses: JamesIves/github-pages-deploy-action@4.1.5
      if: github.event_name != 'pull_request'
      with:
        BRANCH: ${{ env.PROJECT_NAME }}-fpm
        FOLDER: mfi-fpm


