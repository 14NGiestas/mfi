# Contributing to MFI

## Development Setup

### Requirements

- fpm ≥ 0.13.0
- fypp (any version)
- gfortran 12+ (recommended)
- BLAS and LAPACK libraries

For GPU development, additionally:
- CUDA Toolkit (12.x recommended)
- Nix (for reproducible environments)

### Quick Start

```sh
git clone https://github.com/14NGiestas/mfi.git
cd mfi
make              # generates .f90 from .fpp/.fypp templates
fpm test          # runs the test suite
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

### Nix Environments (optional)

For reproducible builds with all dependencies:

```sh
# CPU-only
nix-shell shells/cpu-only.nix --run "fpm test"

# GPU with modern CUDA
nix-shell shells/gpu-modern.nix --run "fpm test --profile cublas"
```

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

The CI pipeline runs automatically only on `main` and PRs targeting `main` to conserve GitHub Actions minutes.

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
| gpu-modern | CPU BLAS/LAPACK (CUDA env) | `release` | Yes |
| gpu-modern-cublas | cuBLAS | `cublas` | Compile only |
| consumer-cpu | CPU consumer project | *(default)* | Build only |
| consumer-cublas | GPU consumer project | `cublas` | Build only |

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
