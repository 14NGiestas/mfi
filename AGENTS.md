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
# Enter dev shell (provides gfortran, fpm, fypp, BLAS, LAPACK)
nix develop             # cpu-only
nix develop .#gpu-modern   # CUDA 12.3
nix develop .#gpu-legacy   # CUDA 11.8

# CPU-only (default)
make
fpm test

# GPU/cuBLAS
make
fpm build --profile cublas
fpm test --profile cublas

# AMD GPU with native HIP and rocBLAS
# Set CPATH and LIBRARY_PATH first. See README.md.
make
fpm build --profile hip
fpm test --profile hip
```

> **Note:** `gpu-modern` CI uses `--profile debug` to avoid a gfortran -O2 optimizer
> bug (fixed in gfortran 15.2.0). See `BUGS.md`.

## Nix Flake

A single `flake.nix` provides all dev shells (replaces old `shells/*.nix`):

| Shell | Command | CUDA |
|-------|---------|------|
| `cpu-only` | `nix develop .#cpu-only` | — |
| `gpu-modern` | `nix develop .#gpu-modern` | 12.3 |
| `gpu-legacy` | `nix develop .#gpu-legacy` | 11.8 |
| default | `nix develop` | — (same as cpu-only) |

- nixpkgs pinned to **24.11** (last version with CUDA 11.8/12.3)
- fpm 0.13.0 via inline overlay (PR #506818 in nixpkgs) — remove once merged
- gfortran, fpm, fypp, pkg-config all provided by the flake
- CI uses `magic-nix-cache-action` for fast cached builds
- Temp make files (`.mfi_*`, `.f77_*, *.tmp`) are gitignored

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
- `hip.fpp` — HIP/rocBLAS C-interop interfaces and constants, plus the `hip_` prefixed versions of the allocation macros
- `src/c/hip_wrap.c` — rocBLAS handle pool, HIP memory helpers, the `MFI_USE_HIP` env var
- `extensions.fpp` — cuBLAS handle lifecycle (`mfi_cublas_handle_ensure`, `mfi_cublas_finalize`), execution mode control (`mfi_force_gpu`, `mfi_force_cpu`), `mfi_cublas_error`, and the rocBLAS equivalents (`mfi_hip_finalize`, `mfi_hip_error`)

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
| `MFI_HIP` | Preprocessor macro | Enables the rocBLAS code at compile time (set by the fpm `hip` feature) |
| `MFI_USE_HIP` | Internal variable + env var | Runtime rocBLAS dispatch flag (read from env var on lazy init) |
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

## HIP and rocBLAS Specifics

- `hip.fpp` and `src/c/hip_wrap.c` mirror `cublas.fpp` and `src/c/cublas_wrap.c`.
- Fypp holds one macro table for the whole run. Both backend files pass through
  fypp. So every HIP macro carries the `hip_` prefix. A second macro named
  `declare` would overwrite the cuBLAS macro.
- rocBLAS has no v2 symbol suffix. The Fortran interfaces bind to `rocblas_sgemm`
  and its siblings in module `f77_blas_hip`.
- rocBLAS enums differ from cuBLAS. `rocblas_operation_none` is 111 and
  `rocblas_operation_transpose` is 112. The constants live in `hip.fpp`.
- The TRSM fill mode needs no inversion here. `rocblas_fill_upper` is 121.
- The wrapper creates one handle per thread and sets the host pointer mode. It
  reads `OMP_NUM_THREADS` for the handle count.
- `rocblas_create_handle` needs a live HIP runtime. `hip_wrap.c` calls
  `hipInit(0)` after a failed create, then it retries once.
- The feature links `librocblas`, `libamdhip64` and `libgomp`.
- `mfi_gemm`, `mfi_gemv`, `mfi_trmm` and `mfi_trsm` check the array shapes in a
  `hip` build. On a mismatch the code calls `error stop` with a message. These
  routines have no status argument, so a caller cannot read an error code. Do
  not turn the check into a returned code.
- The check covers GPU calls and the CPU fallback of a `hip` build.
- Set `CPATH` for the C wrapper and `LIBRARY_PATH` for the link step. Both need
  the rocBLAS and HIP directories. Set `LD_LIBRARY_PATH` for the run.
- A path dependency can lose its `features` entry on a second build with a warm
  fpm cache. Then the app silently uses the CPU stubs. Run a clean build when you
  change a feature list. This hits the cuBLAS consumer projects as well.

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

# AMD GPU with native rocBLAS (set CPATH and LIBRARY_PATH before the build)
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-fpm", features = ["hip"] }
```

**Note:** When using fpm >= 0.13.0, the `cublas` and `cudart` linking requirements are automatically propagated to the consuming project via the `cublas` feature. No explicit `[build] link = [...]` is required in the consumer's `fpm.toml`.

## Testing Notes

- LAPACK tests have pre-existing failures unrelated to cuBLAS work (`cunmrq`, `sormrq`, `heevx` segfault)
- GPU testing available via `gpu_test.ipynb` (Colab: Tesla T4, CUDA 12.8)
- fpm ≥0.13.0 required for `[profiles]` and `[features]` support
- CI uses `MFI_TEST_ELEMENTS=50000` and `MFI_TEST_SAMPLES=1` for fast runs

### Pending test coverage

These wrappers exist but lack tests because the factorization routines (`geqr2`, `gerq2`)
are not yet implemented — they need `tau` from prior factorization:

| Routine | Needs | Status |
|---------|-------|--------|
| `org2r`, `orgr2` | `sgeqr2`/`dgeqr2` | Macro ready, not in `_COLLECT` |
| `ung2r`, `ungr2` | `cgeqr2`/`zgeqr2` | Macro ready, not in `_COLLECT` |
| `orm2r`, `ormr2` | `sgeqr2`/`dgeqr2` | Macro exists, not in `_COLLECT` |
| `unm2r`, `unmr2` | `cgeqr2`/`zgeqr2` | Macro exists, not in `_COLLECT` |

Macros are in `test/lapack/macros/` and can be added to `test/lapack.fpp` once the factorization routines are available.

## fpm.toml Configuration

```toml
[preprocess.cpp]
macros = ["MFI_EXTENSIONS", "MFI_LINK_EXTERNAL"]

[features]
cublas.build.link = ["blas", "lapack", "cublas", "cudart"]
cublas.preprocess.cpp.macros = ["MFI_CUBLAS"]
hip.build.link = ["blas", "lapack", "rocblas", "amdhip64", "gomp"]
hip.preprocess.cpp.macros = ["MFI_HIP"]

[profiles]
cublas = ["cublas"]
hip = ["hip"]
```

- CPU builds use default macros (`MFI_EXTENSIONS`, `MFI_LINK_EXTERNAL`) — no CUDA dependencies.
- `fpm build --profile cublas` activates the `cublas` feature, which adds `MFI_CUBLAS` to the preprocessor and links CUDA libraries.
- Do NOT use `[build] link = [...]` at the root level with `cublas`/`cudart` — those must be gated behind a feature.
