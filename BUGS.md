# Known Issues

## gfortran -O2 optimizer bug with external LAPACK wrappers

**Symptom:** `abs(A_in - A_rf)` returns stale values (the original matrix A) instead of
zeros, even when `A_in` and `A_rf` are bit-identical after LAPACK factorization.
The bug **only** occurs with `-O2` (`--profile release`), not with `-O0` (`--profile debug`).

**Reproducer:**
```fortran
program gfortran_opt_bug
    use f77_lapack, only: spotrf
    use mfi_lapack, only: mfi_spotrf
    use iso_fortran_env, only: real32
    implicit none
    integer, parameter :: wp = real32, N = 3
    real(wp) :: A(N,N), A_rf(N,N), A_in(N,N)
    real(wp) :: diff(N,N)
    integer :: info_rf, info_mfi, info

    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    A_rf = A
    call spotrf('U', N, A_rf, N, info)

    A_in = A
    call mfi_spotrf(A_in, info=info_mfi)

    diff = abs(A_in - A_rf)

    ! With -O0: maxval(diff) = 0.0 ✅
    ! With -O2: maxval(diff) = 2.10263348 ❌ (same as original A(3,3))
    ! Yet printed values of A_rf and A_in are identical!

    if (all(diff < 1e-6_wp)) then
        print '(A)', "PASS"
    else
        print '(A)', "FAIL"
    end if
end program
```

**Root cause:** Likely a register/memory aliasing issue where gfortran -O2 cannot track that
the arrays are modified through external LAPACK calls (which lack `!GCC$ ATTRIBUTES NOALIAS`
or equivalent). The optimizer reuses stale register values instead of reloading from memory.

**Workaround:** Use `--profile debug` (or `-O0`) for testing in CI.

**Affected:** gfortran 13/14 with `-O2` optimization on calls to LAPACK wrappers with
`intent(inout)` array arguments.
