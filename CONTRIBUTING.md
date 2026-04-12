# Contributing to MFI

## Development Setup

### Recommended: Nix Flake (zero config)

```sh
nix develop             # cpu-only shell (gfortran, fpm 0.13, fypp, BLAS, LAPACK)
nix develop .#gpu-modern   # with CUDA 12.3
nix develop .#gpu-legacy   # with CUDA 11.8
```

Requires [Nix](https://nixos.org/download/) with flakes enabled.

### Manual Setup

- fpm ≥ 0.13.0
- fypp (any version)
- gfortran 12+ (recommended)
- BLAS and LAPACK libraries

For GPU development, additionally:
- CUDA Toolkit (12.3 for gpu-modern, 11.8 for gpu-legacy)

### Quick Start

```sh
git clone https://github.com/14NGiestas/mfi.git
cd mfi
nix develop          # enter dev shell (or install deps manually)
make                 # generates .f90 from .fpp/.fypp templates
fpm test             # runs the test suite
```

## Code Generation Architecture

MFI uses [fypp](https://github.com/aradi/fypp) to generate type-generic code.

**The pipeline is:**

```
.fpp / .fypp templates → fypp → csplit (for multi-module files) → .f90 sources → fpm compile
```

**Critical rules:**
- **NEVER edit `.f90` files directly** — they are generated artifacts
- Always modify `.fpp`/`.fypp` templates, then run `make clean && make`
- `.f90` files are gitignored but committed to deployment branches by CI

### File Structure

| Pattern | Purpose |
|---------|---------|
| `*.fpp` / `*.fypp` | Templates (edit these) |
| `src/mfi/blas/*.fypp` | Modern BLAS wrapper implementations |
| `src/f77/blas/*.fypp` | F77 interface declarations |
| `extensions.fpp` | cuBLAS handle lifecycle, execution mode control |
| `cublas.fpp` | CUDA/cuBLAS v2 C-interop interfaces |
| `common.fpp` | Core fypp macros (type prefixes, optional args, interface generators) |

## Build & Test

```sh
# CPU-only (default)
make
fpm test

# GPU/cuBLAS
make
fpm build --profile cublas
fpm test --profile cublas
```

### Nix Environments

The project uses a single `flake.nix` that provides all dev shells:

```sh
nix develop              # cpu-only
nix develop .#gpu-modern # with CUDA 12.3
nix develop .#gpu-legacy # with CUDA 11.8
```

The flake pins nixpkgs to 24.11 (last version with CUDA 11.8/12.3) and provides
fpm 0.13.0 via an overlay until PR #506818 is merged.

## Branch Model

| Branch | Purpose | Deploy Target |
|--------|---------|---------------|
| `main` | Primary development (CPU + cuBLAS via features) | `mfi-fpm` (via CI) |
| `impl/cublas` | GPU/experimental staging | `mfi-cublas` (via CI) |
| `mfi-fpm` | CPU-only deploy artifact (`.f90` + `.toml` only) | — |
| `mfi-cublas` | GPU deploy artifact (`.f90` + `.toml` only) | — |

**CI triggers:**
- Push to `main` → full test matrix → deploy to `mfi-fpm`
- Push to `impl/cublas` → full test matrix → deploy to `mfi-cublas`
- PR to `main` → full test matrix (no deploy)
- Other branches → manual dispatch only

**To deploy:** commit changes and push to the corresponding branch. CI handles the rest.

## Continuous Integration

The CI pipeline uses Nix flakes with `magic-nix-cache-action` for fast, cached builds.

### Manual CI Trigger

For feature branches:
1. Go to the **Actions** tab on GitHub
2. Select the **fpm-deployment** workflow
3. Click **Run workflow** and select your branch
4. The pipeline tests: CPU-only, GPU-modern (with tests), GPU-cuBLAS (compile-only), and consumer projects

### Test Matrix

| Job | Configuration | Profile | Test |
|-----|--------------|---------|------|
| cpu-only | CPU BLAS/LAPACK | *(default)* | Yes |
| gpu-modern | CPU BLAS/LAPACK (CUDA env) | `debug` | Yes |
| gpu-modern-cublas | cuBLAS | `cublas` | Compile only |
| consumer-cpu | CPU consumer project | *(default)* | Build only |
| consumer-cublas | GPU consumer project | `cublas` | Build only |

> **Note:** gpu-modern uses `--profile debug` (not `release`) to avoid a gfortran -O2
> optimizer bug that affects LAPACK wrapper tests. Fixed in gfortran 15.2.0; see `BUGS.md`.

## Conventions

### Naming

| Name | Kind | Purpose |
|------|------|---------|
| `MFI_CUBLAS` | Preprocessor macro | Enables cuBLAS code at compile time (set by fpm `cublas` feature) |
| `MFI_USE_CUBLAS` | Fortran variable + env var | Runtime GPU dispatch flag (set by `mfi_execution_init` from env) |
| `MFI_EXTENSIONS` | Preprocessor macro | Enables BLAS extension routines (iamin, iamax, lamch) |
| `MFI_LINK_EXTERNAL` | Preprocessor macro | Links external BLAS extensions |

### Preprocessor Guards

All CUDA-dependent code must be wrapped:

```fortran
#if defined(MFI_CUBLAS)
    ! cuBLAS-specific code
#else
    ! CPU fallback
#endif
```

### Purity

All MFI BLAS wrappers (`mfi_gemm`, `mfi_gemv`, `mfi_trsm`, `mfi_trmm`) and all `bind(c)` CUDA/cuBLAS interfaces are `pure`. Do NOT remove `pure` — dependent projects call these from `pure` contexts.

### cuBLAS v2 Rules

- All `bind(c)` interfaces must be `pure` with `VALUE` on every argument
- Interface bodies at module level inherit from `use iso_c_binding` — do NOT add `import` inside `interface` blocks
- TRSM fill mode is inverted: `CUBLAS_TRSM_FILL_UPPER = 1`, `CUBLAS_TRSM_FILL_LOWER = 0`
- Use cuBLAS v2 (`cublasSgemm_v2`, etc.), not v1 (`cublasSgemm`)

## Testing Notes

- Known pre-existing LAPACK test failures (unrelated to cuBLAS): `cunmrq`, `sorg2r`, `sorgr2`, `cungr2`, `cung2r`, `sormrq`, `heevx` (segfault)
- GPU testing available via `gpu_test.ipynb` (Colab: Tesla T4, CUDA 12.8)
- CI uses `MFI_TEST_ELEMENTS=50000` and `MFI_TEST_SAMPLES=1` for fast runs
- Override locally: `MFI_TEST_ELEMENTS=1000000 ./build/app/app` for benchmarking

## Dependency Usage

Projects consuming MFI:

```toml
# CPU-only (stable)
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-fpm" }

# GPU/cuBLAS (stable)
mfi = { git="https://github.com/14NGiestas/mfi.git", branch="mfi-cublas", features = ["cublas"] }

# Consuming project must link CUDA when using cublas feature
[build]
link = ["cublas", "cudart"]
```
