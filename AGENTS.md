# MFI — Agent Instructions

## Critical Rules

- **NEVER edit `.f90` files directly** — they are generated artifacts from `.fpp`/`.fypp` sources.
- Always modify `.fpp`/`.fypp` macros, then run:
  ```sh
  make clean && make
  ```
- `.f90` files are gitignored (line 44 of `.gitignore`) but committed to deployment branches by CI.

## Build & Test Commands

```sh
# CPU-only (default)
make
fpm test

# GPU/cuBLAS
make
fpm build --profile cublas
fpm test --profile cublas
```

## Branch Model

| Branch | Purpose | Deployment Target |
|--------|---------|-------------------|
| `main` | Primary development (CPU + cuBLAS via features) | `mfi-fpm` (via CI) |
| `impl/cublas` | GPU/experimental staging | `mfi-cublas` (via CI) |
| `mfi-fpm` | CPU-only deploy artifact (`.f90` + `.toml` only) | — |
| `mfi-cublas` | GPU deploy artifact (`.f90` + `.toml` only) | — |

**CI triggers:**
- Push to `main` → full test matrix → deploy to `mfi-fpm`
- Push to `impl/cublas` → full test matrix → deploy to `mfi-cublas`
- PR to `main` → full test matrix (no deploy)
- Other branches → manual dispatch only

**To deploy:** commit changes, push to the corresponding branch. CI handles the rest.

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

## Naming Conventions

| Name | Kind | Purpose |
|------|------|---------|
| `MFI_CUBLAS` | Preprocessor macro | Enables cuBLAS code at compile time (set by fpm `cublas` feature) |
| `MFI_USE_CUBLAS` | Internal variable + env var | Runtime GPU dispatch flag (read from env var on lazy init) |
| `MFI_EXTENSIONS` | Preprocessor macro | Enables BLAS extension routines (iamin, iamax, lamch) |
| `MFI_LINK_EXTERNAL` | Preprocessor macro | Links external BLAS extensions |

## Purity: Why `pure` on GPU Wrappers is Correct

**All MFI BLAS wrappers (`mfi_gemm`, `mfi_gemv`, `mfi_trsm`, `mfi_trmm`) and all `bind(c)` CUDA/cuBLAS interfaces are `pure`. This is intentional and semantically correct.**

Do NOT remove `pure` from these routines. Reasons:

1. **`error stop` is allowed in `pure` procedures** — permitted by Fortran 2008. The fact that a routine may abort on failure does not make it impure.

2. **GPU alloc → compute → dealloc is semantically pure from Fortran's perspective** — The CUDA device memory (allocated via `cudaMalloc`, freed via `cudaFree`) is opaque to the Fortran compiler. No Fortran-visible state is modified. This is exactly the same pattern as local `allocate`/`deallocate` inside a `pure` CPU function.

3. **`bind(c)` + `pure` is valid** — The compiler cannot verify purity of external C code, so it trusts the declaration. That's the whole point: you're asserting to the compiler that the side effects are not observable from Fortran.

4. **Dependent projects need this** — Projects like CheesyHam call these wrappers from `pure` contexts. Removing `pure` breaks their compilation.

The pattern:
```fortran
pure subroutine mfi_sgemm(a, b, c, ...)
    ! allocate GPU memory (opaque to Fortran)
    ! call cuBLAS (external C, compiler trusts purity claim)
    ! copy result back (no Fortran state modification)
    ! free GPU memory (opaque to Fortran)
    ! return — no observable side effects
end subroutine
```
is identical to:
```fortran
pure function foo(x) result(y)
    real, allocatable :: tmp(:)
    allocate(tmp(size(x)))    ! allowed in pure
    tmp = x * 2.0
    y = sum(tmp)
    deallocate(tmp)           ! allowed in pure
end function
```

## cuBLAS v2 Specifics

- **All `bind(c)` interfaces must be `pure`** with `VALUE` on every argument (including `intent(out)` pointers)
- **All MFI BLAS wrappers (`mfi_gemm`, `mfi_gemv`, `mfi_trsm`, `mfi_trmm`) must be `pure`**
- **Interface bodies at module level inherit from `use iso_c_binding`** — do NOT add `import`, `use`, or `import ::` inside `interface` blocks at module scope. Host association handles type visibility.
- cuBLAS stat checks use `call mfi_cublas_error(stat, 'name')` (a pure subroutine wrapper) for consistency with the purity design
- **TRSM fill mode is inverted:** `CUBLAS_TRSM_FILL_UPPER = 1`, `CUBLAS_TRSM_FILL_LOWER = 0` (opposite of standard BLAS enums)
- cuBLAS v1 (`cublasAlloc`/`cublasSgemm`) is deprecated — use v2 (`cudaMalloc`, `cublasCreate_v2`, `cublasSgemm_v2`, etc.)

## Runtime CPU/GPU Switching

- **Zero-config default:** `call mfi_gemm(A, B, C)` — always works, no setup
- **Env var activation:** `MFI_USE_CUBLAS=1 ./app` — lazy init reads env var automatically
- **Manual switch:** `call mfi_force_gpu` / `call mfi_force_cpu` — **always available** (stub no-op when compiled without `cublas`, functional with `cublas` feature)
- **OpenMP safe:** Per-thread cuBLAS handles, pre-allocated from `OMP_NUM_THREADS`
- **No state leaks:** Each `mfi_force_*` resets lazy-init state, so calls are safe to repeat

## Dependency Usage

```toml
# CPU-only (stable)
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-fpm" }

# GPU/cuBLAS (stable)
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-cublas", features = ["cublas"] }
```

The consuming project must also link CUDA libraries when using the `cublas` feature:
```toml
[build]
link = ["cublas", "cudart"]
```

## Testing Notes

- LAPACK tests have pre-existing failures unrelated to cuBLAS work (`cunmrq`, `sorg2r`, `sorgr2`, `cungr2`, `cung2r`, `sormrq`, `heevx` segfault)
- GPU testing available via `gpu_test.ipynb` (Colab: Tesla T4, CUDA 12.8)
- fpm ≥0.13.0 required for `[profiles]` and `[features]` support

## fpm.toml Configuration

```toml
[preprocess.cpp]
macros = ["MFI_EXTENSIONS", "MFI_LINK_EXTERNAL"]

[features]
cublas.build.link = ["blas", "lapack", "cublas", "cudart"]
cublas.preprocess.cpp.macros = ["MFI_CUBLAS"]

[profiles]
cublas = ["cublas"]
```

- CPU builds use default macros (`MFI_EXTENSIONS`, `MFI_LINK_EXTERNAL`) — no CUDA dependencies.
- `fpm build --profile cublas` activates the `cublas` feature, which adds `MFI_CUBLAS` to the preprocessor and links CUDA libraries.
- Do NOT use `[build] link = [...]` at the root level with `cublas`/`cudart` — those must be gated behind a feature.
