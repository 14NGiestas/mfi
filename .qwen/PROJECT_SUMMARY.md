# Project Summary

## Overall Goal
Restructure the MFI (Modern Fortran Interfaces) project test suite to make BLAS and LAPACK tests follow consistent patterns, moving from consolidated test files to individual test files for better maintainability.

## Key Knowledge
- Project is located at `/home/ian/gits/mfi/` with Fortran source code in `src/mfi/` and legacy F77 code in `src/f77/`
- Uses fpm (Fortran Package Manager) with `fpm.toml` configuration and Makefile for build process
- Built with gfortran compiler, with fypp (Fortran preprocessor) for generating interfaces
- BLAS tests were originally in a single consolidated file `test/blas.fpp`, while LAPACK tests were already in individual files
- Uses conditional compilation with `-DMFI_EXTENSIONS` flag for BLAS extension routines (iamin/iamax)
- Both BLAS and LAPACK tests now use consistent structure with `test_run` macro generating timed test blocks
- Shared `test/assert.inc` file provides consistent assert subroutine with optional info parameter
- Tests are built with `make` and run with `fpm test`, with extensions enabled via `FYPPFLAGS=-DMFI_EXTENSIONS`

## Recent Actions
- [DONE] Converted single BLAS test file (`test/blas.fpp`) to 16 individual test files: asum.fpp, axpy.fpp, copy.fpp, dot.fpp, dotc.fpp, dotu.fpp, gemm.fpp, gemv.fpp, iamax.fpp, iamin.fpp, lamch.fpp, nrm2.fpp, rot.fpp, rotg.fpp, scal.fpp, swap.fpp
- [DONE] Moved BLAS test macros to `test/blas/macros/` directory for better organization
- [DONE] Added conditional compilation to extension tests (iamin, iamax) to only run when MFI_EXTENSIONS is defined
- [DONE] Created shared `test/assert.inc` for common assert subroutine
- [DONE] Made LAPACK tests consistent with BLAS tests by updating all LAPACK test files to use the same `test_run` macro structure
- [DONE] Enhanced assert subroutine to support optional info parameter for better error reporting
- [DONE] Verified all changes maintain compatibility with both Makefile and fpm build systems

## Current Plan
- [DONE] Restructure BLAS tests from consolidated to individual files
- [DONE] Make LAPACK test structure consistent with BLAS test structure  
- [DONE] Ensure build system compatibility with both approaches
- [TODO] Address pre-existing LAPACK hetrf test failure (unrelated to restructuring changes)
- [TODO] Consider if additional test consistency improvements are needed

---

## Summary Metadata
**Update time**: 2025-12-10T21:57:32.317Z 
