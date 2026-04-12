# MFI

## Modern Fortran interfaces to BLAS and LAPACK

MFI provides generic, type-agnostic wrappers around BLAS and LAPACK routines.
Instead of writing type-specific calls with dozens of arguments, you write one
call that works for `real32`, `real64`, `complex(real32)`, and `complex(real64)`.

### Example: $C = A \cdot B$

```fortran
program main
    use mfi_blas, only: mfi_gemm
    implicit none
    real :: A(4,4), B(4,4), C(4,4)
    ! ... fill A and B ...
    call mfi_gemm(A, B, C)   ! That's it. No leading dims, no m/n/k, no alpha/beta.
end program
```

---

## Quick Start

### Recommended: Nix Flake (zero config)

```sh
git clone https://github.com/14NGiestas/mfi.git
cd mfi
nix develop          # cpu-only shell with gfortran, fpm, fypp, BLAS, LAPACK
nix develop .#gpu-modern   # with CUDA 12.3
nix develop .#gpu-legacy   # with CUDA 11.8
make              # generates .f90 from .fpp/.fypp templates
fpm test          # runs the test suite
```

Requires [Nix](https://nixos.org/download/) with flakes enabled.

### Manual setup

| Tool | Minimum version |
|------|-----------------|
| fpm | ≥ 0.13.0 |
| fypp | any |
| Fortran compiler | gfortran 12+ (recommended) |

```sh
pip install fypp
```

Install BLAS and LAPACK from your package manager:

| Distro | Package |
|--------|---------|
| Arch | `openblas-lapack-static` (AUR) |
| Ubuntu/Debian | `libblas-dev liblapack-dev` |
| Fedora | `openblas-devel lapack-devel` |

### Build & Test

```sh
git clone https://github.com/14NGiestas/mfi.git
cd mfi
make              # generates .f90 from .fpp/.fypp templates
fpm test          # runs the test suite
```

---

## Using MFI as a Dependency

Add to your project's `fpm.toml`:

```toml
# CPU-only (stable)
[dependencies]
mfi = { git = "https://github.com/14NGiestas/mfi.git", branch = "mfi-fpm" }
```

That's all — fpm handles the rest. No `make` needed in your own project.

---

## GPU Acceleration with cuBLAS

MFI can transparently dispatch BLAS calls to cuBLAS when compiled with the
`cublas` feature. The same `mfi_gemm`, `mfi_gemv`, etc. calls run on the GPU
without code changes.

Try it in your browser:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/14NGiestas/mfi/blob/main/gpu_test.ipynb)

### Local build with cuBLAS

```sh
make
fpm build --profile cublas
fpm test --profile cublas
```

### Runtime CPU / GPU switching

MFI uses lazy initialization — no setup code is needed. When compiled with the
`cublas` feature, GPU dispatching is controlled entirely by the
`MFI_USE_CUBLAS` environment variable:

```sh
# CPU (default)
./build/app/app

# GPU
MFI_USE_CUBLAS=1 ./build/app/app
```

The same `call mfi_gemm(A, B, C)` runs on CPU or GPU without any code changes.

For OpenMP-parallel programs, also set `OMP_NUM_THREADS` to pre-allocate
per-thread cuBLAS handles:

```sh
MFI_USE_CUBLAS=1 OMP_NUM_THREADS=8 ./build/app/app
```

#### Manual CPU/GPU switching (advanced)

If you need fine-grained control within a single program (e.g., run most
computations on GPU but force a specific call to CPU), use
`mfi_force_gpu` / `mfi_force_cpu`:

```fortran
call mfi_gemm(A, B, C)       ! CPU (default)

call mfi_force_gpu
call mfi_gemm(D, E, F)       ! GPU
call mfi_force_cpu

call mfi_gemm(G, H, I)       ! CPU again
```

> **Note:** When compiled **without** the `cublas` feature, `mfi_force_gpu` and
> `mfi_force_cpu` are no-op stubs — your code compiles and runs normally on CPU
> without any `#ifdef` changes. Simply recompile with `--profile cublas` to
> activate GPU acceleration.

#### Clean shutdown (optional)

Call `mfi_cublas_finalize()` at program end to release GPU resources.
The OS cleans up on exit anyway.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `CUBLAS_STATUS_NOT_INITIALIZED` | cuBLAS handle not created. Set `MFI_USE_CUBLAS=1` or call `mfi_force_gpu` before the first BLAS call. |
| `cuda_runtime.h not found` | CUDA Toolkit is not installed or not in your include path. See [gpu_test.ipynb](gpu_test.ipynb) for a working Colab setup. |
| `i?amin` symbols missing | Your BLAS provider lacks extensions. Use the default profile (without `MFI_LINK_EXTERNAL`) or switch to OpenBLAS. |
| Tests fail on CPU build | Known pre-existing failures: `cunmrq`, `sorg2r`, `sorgr2`, `cungr2`, `cung2r`, `sormrq`, `heevx` (segfault). |

---

## Interface Levels

MFI exposes four interface levels for BLAS, from bare-metal to fully modern:

| Level | Example | Arguments |
|-------|---------|-----------|
| Raw F77 | `call cgemm('N','N', N, N, N, alpha, A, N, B, N, beta, C, N)` | 13 |
| Improved F77 | `call f77_gemm('N','N', N, N, N, alpha, A, N, B, N, beta, C, N)` | 13 (no `c`/`d`/`s`/`z` prefix) |
| MFI typed | `call mfi_sgemm(A, B, C)` | 3 (type-specific) |
| MFI generic | `call mfi_gemm(A, B, C)` | 3 (type-agnostic) |

For full API documentation, see the [generated reference](https://14ngiestas.github.io/mfi/).

---

## Supported Routines

### BLAS

#### Level 1

<details>
<summary>Click to expand</summary>

| Status | Name   | Description |
|--------|--------|-------------|
| :+1: | asum   | Sum of vector magnitudes |
| :+1: | axpy   | Scalar-vector product |
| :+1: | copy   | Copy vector |
| :+1: | dot    | Dot product |
| :+1: | dotc   | Dot product conjugated |
| :+1: | dotu   | Dot product unconjugated |
| f77 | sdsdot | Extended precision inner product |
| f77 | dsdot  | Extended precision inner product with double result |
| :+1: | nrm2   | Vector 2-norm (Euclidean norm) |
| :+1: | rot    | Plane rotation |
| :+1: | rotg   | Generate Givens rotation |
| :+1: | rotm   | Modified Givens rotation |
| :+1: | rotmg  | Generate modified Givens rotation |
| :+1: | scal   | Vector-scalar product |
| :+1: | swap   | Vector-vector swap |

</details>

#### Level 1 — Extensions

<details>
<summary>Click to expand</summary>

| Status | Name  | Description |
|--------|-------|-------------|
| :+1: | iamax | Index of maximum absolute value element |
| :+1: | iamin | Index of minimum absolute value element |
| :+1: | lamch | Machine precision parameters |

</details>

#### Level 2

<details>
<summary>Click to expand</summary>

| Status | Name | Description |
|--------|------|-------------|
| :+1: | gbmv | Matrix-vector product (general band) |
| :+1: | gemv | Matrix-vector product (general) |
| :+1: | ger  | Rank-1 update (general) |
| :+1: | gerc | Rank-1 update (general, conjugated) |
| :+1: | geru | Rank-1 update (general, unconjugated) |
| :+1: | hbmv | Matrix-vector product (Hermitian band) |
| :+1: | hemv | Matrix-vector product (Hermitian) |
| :+1: | her  | Rank-1 update (Hermitian) |
| :+1: | her2 | Rank-2 update (Hermitian) |
| :+1: | hpmv | Matrix-vector product (Hermitian packed) |
| :+1: | hpr  | Rank-1 update (Hermitian packed) |
| :+1: | hpr2 | Rank-2 update (Hermitian packed) |
| :+1: | sbmv | Matrix-vector product (symmetric band) |
| :+1: | spmv | Matrix-vector product (symmetric packed) |
| :+1: | spr  | Rank-1 update (symmetric packed) |
| :+1: | spr2 | Rank-2 update (symmetric packed) |
| :+1: | symv | Matrix-vector product (symmetric) |
| :+1: | syr  | Rank-1 update (symmetric) |
| :+1: | syr2 | Rank-2 update (symmetric) |
| :+1: | tbmv | Matrix-vector product (triangular band) |
| :+1: | tbsv | Solve (triangular band) |
| :+1: | tpmv | Matrix-vector product (triangular packed) |
| :+1: | tpsv | Solve (triangular packed) |
| :+1: | trmv | Matrix-vector product (triangular) |
| :+1: | trsv | Solve (triangular) |

</details>

#### Level 3

<details>
<summary>Click to expand</summary>

| Status | GPU | Name  | Description |
|--------|-----|-------|-------------|
| :+1: | :white_check_mark: | gemm  | General matrix-matrix product |
| :+1: | :white_check_mark: | hemm  | Hermitian × general matrix product |
| :+1: | | herk  | Hermitian rank-k update |
| :+1: | | her2k | Hermitian rank-2k update |
| :+1: | :white_check_mark: | symm  | Symmetric × general matrix product |
| :+1: | | syrk  | Symmetric rank-k update |
| :+1: | | syr2k | Symmetric rank-2k update |
| :+1: | :white_check_mark: | trmm  | Triangular × general matrix product |
| :+1: | :white_check_mark: | trsm  | Solve with triangular matrix |

</details>

### LAPACK

LAPACK coverage is growing — routines are implemented as needed.

#### Factorization and Solve

<details>
<summary>Click to expand</summary>

| Status | Name  | Description |
|--------|-------|-------------|
| :+1: | geqrf | QR factorization |
| :+1: | gerqf | RQ factorization |
| :+1: | getrf | LU factorization |
| :+1: | getri | Matrix inverse (from LU) |
| :+1: | getrs | Solve with LU-factored matrix |
| :+1: | hetrf | Bunch-Kaufman factorization (Hermitian) |
| :+1: | pocon | Condition number estimate (Cholesky) |
| :+1: | potrf | Cholesky factorization |
| :+1: | potri | Matrix inverse (from Cholesky) |
| :+1: | potrs | Solve with Cholesky-factored matrix |
| :+1: | sytrf | Bunch-Kaufman factorization (symmetric) |
| :+1: | trtrs | Solve with triangular matrix |

</details>

#### Orthogonal / Unitary Factors

<details>
<summary>Click to expand</summary>

| Status | Name  | Description |
|--------|-------|-------------|
| :+1: | orgqr | Generate Q from QR (real) |
| :+1: | orgrq | Generate Q from RQ (real) |
| :+1: | ormqr | Multiply by Q from QR (real) |
| f77 | ormrq | Multiply by Q from RQ (real) |
| :+1: | org2r | Generate Q from QR2 (real) |
| :+1: | orm2r | Multiply by Q from QR2 (real) |
| :+1: | orgr2 | Generate Q from RQ2 (real) |
| :+1: | ormr2 | Multiply by Q from RQ2 (real) |
| :+1: | ungqr | Generate Q from QR (complex) |
| :+1: | ungrq | Generate Q from RQ (complex) |
| :+1: | unmqr | Multiply by Q from QR (complex) |
| f77 | unmrq | Multiply by Q from RQ (complex) |
| :+1: | ung2r | Generate Q from QR2 (complex) |
| :+1: | unm2r | Multiply by Q from QR2 (complex) |
| :+1: | ungr2 | Generate Q from RQ2 (complex) |
| :+1: | unmr2 | Multiply by Q from RQ2 (complex) |

</details>

#### Eigenvalues and SVD

<details>
<summary>Click to expand</summary>

| Status | Name  | Description |
|--------|-------|-------------|
| :+1: | gesvd | Singular value decomposition |
| :+1: | heevd | Hermitian eigenvalues (divide & conquer) |
| :+1: | hegvd | Generalized Hermitian eigenproblem (divide & conquer) |
| :+1: | heevr | Hermitian eigenvalues (relatively robust) |
| f77 | heevx | Hermitian eigenvalues (expert) |

</details>

#### Least Squares

<details>
<summary>Click to expand</summary>

| Status | Name   | Description |
|--------|--------|-------------|
| f77 | gels   | Least squares (QR/LQ) |
| f77 | gelst  | Least squares (QR/LQ, T matrix) |
| f77 | gelss  | Least squares (SVD, QR iteration) |
| f77 | gelsd  | Least squares (SVD, divide & conquer) |
| f77 | gelsy  | Least squares (complete orthogonal) |
| f77 | getsls | Least squares (tall-skinny QR/LQ) |
| f77 | gglse  | Equality-constrained least squares |
| f77 | ggglm  | Gauss-Markov linear model |

</details>

#### Auxiliary

| Name      | Types | Description |
|-----------|-------|-------------|
| mfi_lartg | s, d, c, z | Generate plane rotation |

---

## Continuous Integration

CI uses [Nix flakes](flake.nix) with `magic-nix-cache-action` for fast, reproducible builds.

| Event | Behavior |
|-------|----------|
| Push to `main` | Full test matrix + deploy to `mfi-fpm` |
| Push to `impl/cublas` | Full test matrix + deploy to `mfi-cublas` |
| PR to `main` | Full test matrix |
| Manual dispatch | Full test matrix |
