# Known Issues

## gfortran -O2 optimizer bug with external LAPACK wrappers ~~(fixed in gfortran 15)~~

**Status:** Fixed in gfortran 15.2.0. No longer affects CI.

**Historical record:** With gfortran 13/14 and `-O2`, `abs(A_in - A_rf)` returned stale values
(the original matrix A) instead of zeros after LAPACK factorization. The optimizer at `-O2`
failed to reload array values from memory after external subroutine calls, reusing stale
register contents instead.

**Workaround (for older gfortran):** Use `--profile debug` (or `-O0`) for testing in CI.

**Affected:** gfortran 13/14 with `-O2` optimization. Not reproducible on gfortran 15.2.0.
