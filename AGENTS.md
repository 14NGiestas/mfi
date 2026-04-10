# MFI — Agent Instructions

## Critical Rules

- **NEVER edit `.f90` files directly** — they are generated artifacts from `.fpp`/`.fypp` sources.
- Always modify `.fpp`/`.fypp` macros, then run:
  ```sh
  make clean && make FYPPFLAGS="-DMFI_EXTENSIONS -DMFI_USE_CUBLAS"
  ```
- `.f90` files are gitignored (line 44 of `.gitignore`) but committed to deployment branches by CI.

## Build & Test Commands

```sh
# CPU-only (default)
make FYPPFLAGS=-DMFI_EXTENSIONS
fpm test

# With BLAS extensions
make FYPPFLAGS="-DMFI_EXTENSIONS -DMFI_LINK_EXTERNAL"
fpm test

# GPU/cuBLAS
make FYPPFLAGS="-DMFI_EXTENSIONS -DMFI_USE_CUBLAS"
MFI_USE_CUBLAS=1 fpm test --profile cublas
```

## Branch Model

| Branch | Purpose | Deployment Target |
|--------|---------|-------------------|
| `main` | CPU-only, rolling release | `mfi-fpm` (via CI) |
| `impl/cublas` | GPU/experimental development | `mfi-cublas` (via CI) |
| `mfi-fpm` | CPU-only deploy artifact (`.f90` + `.toml` only) | — |
| `mfi-cublas` | GPU deploy artifact (`.f90` + `.toml` only) | — |

**CI triggers:**
- Push to `main` → `test-and-deploy-cpu` job → deploys to `mfi-fpm`
- Push to `impl/cublas` → `test-and-deploy-cublas` job → deploys to `mfi-cublas`

**To deploy cuBLAS changes:** merge working branch into `impl/cublas`, commit regenerated `.f90` files, then push.

## Code Generation Architecture

### Macro Files (edit these)
- `common.fpp` — Core fypp macros: type prefixes (`s`,`d`,`c`,`z`), `@:optional`, `@:defaults`, interface generators
- `cublas.fpp` — CUDA/cuBLAS v2 C-interop interfaces (`pure` + `VALUE` on all `bind(c)` args), `@:allocate`, `@:deallocate`, `@:set_matrix`, `@:get_matrix` macros, cuBLAS constants
- `extensions.fpp` — cuBLAS handle lifecycle (`mfi_cublas_handle_ensure`, `mfi_cublas_finalize`), execution mode control (`mfi_force_gpu`, `mfi_force_cpu`), `mfi_cublas_error`

### Source Macros (edit these)
- `src/mfi/blas/*.fypp` — MFI modern wrapper implementations
- `src/f77/blas/*.fypp` — F77 interface declarations
- `src/mfi/lapack/*.fypp` — LAPACK modern wrappers
- `src/f77/lapack/*.fypp` — LAPACK F77 interfaces

### Generated (do not edit)
- `src/f77/blas.f90`, `src/mfi/blas.f90`, `src/f77/lapack.f90`, `src/mfi/lapack.f90`
- All `test/**/*.f90` files

## cuBLAS v2 Specifics

- **All `bind(c)` interfaces must be `pure`** with `VALUE` on every argument (including `intent(out)` pointers)
- **All MFI BLAS wrappers (`mfi_gemm`, `mfi_gemv`, `mfi_trsm`, `mfi_trmm`) must be `pure`**
- `error stop` is allowed in `pure` subroutines per Fortran standard, but cuBLAS stat checks use `call mfi_cublas_error(stat, 'name')` for consistency
- **TRSM fill mode is inverted:** `CUBLAS_TRSM_FILL_UPPER = 1`, `CUBLAS_TRSM_FILL_LOWER = 0` (opposite of standard BLAS enums)
- cuBLAS v1 (`cublasAlloc`/`cublasSgemm`) is deprecated — use v2 (`cudaMalloc`, `cublasCreate_v2`, `cublasSgemm_v2`, etc.)

## Dependency Usage

Projects consuming MFI as a dependency:
```toml
# CPU-only
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-fpm" }

# GPU/cuBLAS
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-cublas" }
```

## Testing Notes

- LAPACK tests have pre-existing failures unrelated to cuBLAS work (cunmrq, sorg2r, sorgr2, cungr2, cung2r, sormrq, heevx segfault)
- GPU testing available via `gpu_test.ipynb` (Colab: Tesla T4, CUDA 12.8)
- fpm ≥0.13.0 required for `[profiles]` and `[features]` support (v0.10.0 in CI)

## fpm.toml Quirks

- `link-flag` in root `[build]` table is rejected by fpm 0.13.0+ — use `[profiles]` and `[features]` instead
- Current config:
  ```toml
  [build]
  link = ["blas","lapack","cublas","cudart"]
  
  [features]
  cublas.flags = "-I/usr/local/cuda/include -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64"
  
  [profiles]
  cublas = ["cublas"]
  ```
