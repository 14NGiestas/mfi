program main
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_slamch 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_lamch against slamch", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dlamch 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_lamch against dlamch", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_scopy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_copy against scopy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dcopy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_copy against dcopy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ccopy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_copy against ccopy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zcopy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_copy against zcopy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against sswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against dswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against cswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against zswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_saxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against saxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_daxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against daxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_caxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against caxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zaxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against zaxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against sgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against dgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against cgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against zgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against sgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against dgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against cgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against zgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_asum against sasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_asum against dasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_scasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_asum against scasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dzasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_asum against dzasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_snrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 against snrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dnrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 against dnrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_scnrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 against scnrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dznrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 against dznrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_srot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against srot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_drot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against drot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_crot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against crot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against zrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_csrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against csrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against zdrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_scal against sscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_scal against dscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_scal against cscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_scal against zscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_csscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_scal against csscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_scal against zdscal", t2-t1
end block

block
real :: t1, t2
call cpu_time(t1)
 call test_isamin 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamin against isamin", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_idamin 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamin against idamin", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_icamin 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamin against icamin", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_izamin 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamin against izamin", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_isamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against isamax", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_idamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against idamax", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_icamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against icamax", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_izamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against izamax", t2-t1
end block

contains

subroutine test_slamch
    use f77_blas, only: slamch
    use mfi_blas, only: mfi_lamch

    integer, parameter :: wp = REAL32
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
    real(wp) :: a, b 
    integer :: i
    
    do i=1,size(options)
        a = slamch(options(i))
        b = mfi_lamch(options(i),1.0_wp)
        call assert(a == b, "different results for option "//options(i))
    end do

end subroutine
subroutine test_dlamch
    use f77_blas, only: dlamch
    use mfi_blas, only: mfi_lamch

    integer, parameter :: wp = REAL64
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
    real(wp) :: a, b 
    integer :: i
    
    do i=1,size(options)
        a = dlamch(options(i))
        b = mfi_lamch(options(i),1.0_wp)
        call assert(a == b, "different results for option "//options(i))
    end do

end subroutine
subroutine test_scopy
    use f77_blas, only: scopy, f77_copy
    use mfi_blas, only: mfi_copy, mfi_scopy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    real(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(X)
    Y = 0.0_wp

    x_in = x
    y_in = y
    ! The test is always against the original
    call scopy(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_copy(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_scopy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_dcopy
    use f77_blas, only: dcopy, f77_copy
    use mfi_blas, only: mfi_copy, mfi_dcopy

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    real(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(X)
    Y = 0.0_wp

    x_in = x
    y_in = y
    ! The test is always against the original
    call dcopy(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_copy(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_dcopy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_ccopy
    use f77_blas, only: ccopy, f77_copy
    use mfi_blas, only: mfi_copy, mfi_ccopy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    y%re = 0.0_wp
    y%im = 0.0_wp

    x_in = x
    y_in = y
    ! The test is always against the original
    call ccopy(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_copy(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_ccopy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_zcopy
    use f77_blas, only: zcopy, f77_copy
    use mfi_blas, only: mfi_copy, mfi_zcopy

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    y%re = 0.0_wp
    y%im = 0.0_wp

    x_in = x
    y_in = y
    ! The test is always against the original
    call zcopy(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_copy(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_zcopy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_sswap
    use f77_blas, only: sswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_sswap

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    real(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(X)
    call random_number(Y)

    x_in = x
    y_in = y
    ! The test is always against the original
    call sswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_sswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_dswap
    use f77_blas, only: dswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_dswap

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    real(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(X)
    call random_number(Y)

    x_in = x
    y_in = y
    ! The test is always against the original
    call dswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_dswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_cswap
    use f77_blas, only: cswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_cswap

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    x_in = x
    y_in = y
    ! The test is always against the original
    call cswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_cswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_zswap
    use f77_blas, only: zswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_zswap

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(wp) :: rnd(N)

    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    x_in = x
    y_in = y
    ! The test is always against the original
    call zswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_zswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_saxpy
    use f77_blas, only: saxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_saxpy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    real(wp) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: alpha

    call random_number(X)
    call random_number(Y)
    call random_number(alpha)

    x_in = X
    y_in = Y
    call saxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_saxpy(x_in,y_in,alpha)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_daxpy
    use f77_blas, only: daxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_daxpy

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    real(wp) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: alpha

    call random_number(X)
    call random_number(Y)
    call random_number(alpha)

    x_in = X
    y_in = Y
    call daxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_daxpy(x_in,y_in,alpha)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_caxpy
    use f77_blas, only: caxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_caxpy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    complex(wp) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    complex(wp) :: alpha

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

    x_in = X
    y_in = Y
    call caxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_caxpy(x_in,y_in,alpha)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_zaxpy
    use f77_blas, only: zaxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_zaxpy

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    complex(wp) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    complex(wp) :: alpha

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

    x_in = X
    y_in = Y
    call zaxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_zaxpy(x_in,y_in,alpha)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_sgemv
    use f77_blas, only: sgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_sgemv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    real(wp) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

    call random_number(M)
    call random_number(X)
    call random_number(Y)
    call random_number(alpha)
    call random_number(beta)


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call sgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_sgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_dgemv
    use f77_blas, only: dgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_dgemv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    real(wp) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

    call random_number(M)
    call random_number(X)
    call random_number(Y)
    call random_number(alpha)
    call random_number(beta)


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call dgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_dgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_cgemv
    use f77_blas, only: cgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_cgemv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(wp) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    complex(wp) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im)
end block
block
    real(wp) :: re(N)
    real(wp) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(wp) :: re(N)
    real(wp) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im)
end block


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call cgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_cgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_zgemv
    use f77_blas, only: zgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_zgemv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(wp) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    complex(wp) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im)
end block
block
    real(wp) :: re(N)
    real(wp) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(wp) :: re(N)
    real(wp) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im)
end block


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call zgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_zgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_sgemm
    use f77_blas, only: sgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_sgemm

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    real(wp) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    real(wp) :: alpha, beta
    character :: transa, transb
    integer :: i, j

    call random_number(A)
    call random_number(B)
    call random_number(C)
    call random_number(alpha)
    call random_number(beta)


    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call sgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_sgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_dgemm
    use f77_blas, only: dgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_dgemm

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    real(wp) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    real(wp) :: alpha, beta
    character :: transa, transb
    integer :: i, j

    call random_number(A)
    call random_number(B)
    call random_number(C)
    call random_number(alpha)
    call random_number(beta)


    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call dgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_dgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_cgemm
    use f77_blas, only: cgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_cgemm

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    complex(wp) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    complex(wp) :: alpha, beta
    character :: transa, transb
    integer :: i, j

block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im)
end block
block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im)
end block
block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    C = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im)
end block


    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call cgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_cgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_zgemm
    use f77_blas, only: zgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_zgemm

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    complex(wp) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    complex(wp) :: alpha, beta
    character :: transa, transb
    integer :: i, j

block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im)
end block
block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im)
end block
block
    real(wp) :: re(N,N)
    real(wp) :: im(N,N)
    call random_number(im)
    call random_number(re)
    C = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(wp) :: re
    real(wp) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im)
end block


    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call zgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_zgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_sasum
    use f77_blas, only: sasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_sasum

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    real(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = sasum(N, array, 1)
    res(2) = mfi_sasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = sasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_sasum(array)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dasum
    use f77_blas, only: dasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_dasum

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    real(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dasum(N, array, 1)
    res(2) = mfi_dasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = dasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_dasum(array)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_scasum
    use f77_blas, only: scasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_scasum

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    complex(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = scasum(N, array, 1)
    res(2) = mfi_scasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = scasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_scasum(array)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dzasum
    use f77_blas, only: dzasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_dzasum

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    complex(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dzasum(N, array, 1)
    res(2) = mfi_dzasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = dzasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_dzasum(array)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_snrm2
    use f77_blas, only: snrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_snrm2

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    real(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = snrm2(N, array, 1)
    res(2) = mfi_snrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = snrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_snrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dnrm2
    use f77_blas, only: dnrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_dnrm2

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    real(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dnrm2(N, array, 1)
    res(2) = mfi_dnrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = dnrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_dnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_scnrm2
    use f77_blas, only: scnrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_scnrm2

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    complex(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = scnrm2(N, array, 1)
    res(2) = mfi_scnrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = scnrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_scnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dznrm2
    use f77_blas, only: dznrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_dznrm2

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    complex(wp) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dznrm2(N, array, 1)
    res(2) = mfi_dznrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = dznrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_dznrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_srot
    use f77_blas, only: srot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_srot

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: angle
    real(wp) :: c
    real(wp) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(X)
    call random_number(Y)

    c = cos(angle)

    s = i * sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call srot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_srot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_drot
    use f77_blas, only: drot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_drot

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: angle
    real(wp) :: c
    real(wp) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(X)
    call random_number(Y)

    c = cos(angle)

    s = i * sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call drot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_drot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_crot
    use f77_blas, only: crot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_crot

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: angle
    real(wp) :: c
    complex(wp) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    c = cos(angle)

    s = i * sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call crot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_crot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_zrot
    use f77_blas, only: zrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_zrot

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: angle
    real(wp) :: c
    complex(wp) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    c = cos(angle)

    s = i * sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call zrot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_zrot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_csrot
    use f77_blas, only: csrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_csrot

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: angle
    real(wp) :: c
    real(wp) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    c = cos(angle)

    s = sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call csrot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_csrot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_zdrot
    use f77_blas, only: zdrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_zdrot

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(wp) :: angle
    real(wp) :: c
    real(wp) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    c = cos(angle)

    s = sin(angle)

    x_in = x
    y_in = y
    ! The test is always against the original
    call zdrot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_zdrot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_sscal
    use f77_blas, only: sscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_sscal

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    real(wp) :: x(N),    &
                x_in(N), &
                x_rf(N)
    real(wp) :: alpha

    call random_number(X)
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call sscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_sscal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_dscal
    use f77_blas, only: dscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_dscal

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    real(wp) :: x(N),    &
                x_in(N), &
                x_rf(N)
    real(wp) :: alpha

    call random_number(X)
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call dscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_dscal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_cscal
    use f77_blas, only: cscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_cscal

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    complex(wp) :: x(N),    &
                x_in(N), &
                x_rf(N)
    complex(wp) :: alpha

    call random_number(rnd_vector)
    x%re = rnd_vector
    call random_number(rnd_vector)
    x%im = rnd_vector
    call random_number(rnd)
    alpha%re = rnd
    call random_number(rnd)
    alpha%im = rnd

    ! The test is always against the original
    x_in = x
    call cscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_cscal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_zscal
    use f77_blas, only: zscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_zscal

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    complex(wp) :: x(N),    &
                x_in(N), &
                x_rf(N)
    complex(wp) :: alpha

    call random_number(rnd_vector)
    x%re = rnd_vector
    call random_number(rnd_vector)
    x%im = rnd_vector
    call random_number(rnd)
    alpha%re = rnd
    call random_number(rnd)
    alpha%im = rnd

    ! The test is always against the original
    x_in = x
    call zscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_zscal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_csscal
    use f77_blas, only: csscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_csscal

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    complex(wp) :: x(N),    &
                x_in(N), &
                x_rf(N)
    real(wp) :: alpha

    call random_number(rnd_vector)
    x%re = rnd_vector
    call random_number(rnd_vector)
    x%im = rnd_vector
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call csscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_csscal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_zdscal
    use f77_blas, only: zdscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_zdscal

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(wp) :: rnd_vector(N), rnd
    complex(wp) :: x(N),    &
                x_in(N), &
                x_rf(N)
    real(wp) :: alpha

    call random_number(rnd_vector)
    x%re = rnd_vector
    call random_number(rnd_vector)
    x%im = rnd_vector
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call zdscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_zdscal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(x_in, alpha)
    call assert(all(x_in == x_rf), "different results")

end subroutine

subroutine test_isamin
    use f77_blas, only: isamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_isamin

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = isamin(N, array, 1)
    res(2) = mfi_isamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = isamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_isamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_idamin
    use f77_blas, only: idamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_idamin

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = idamin(N, array, 1)
    res(2) = mfi_idamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = idamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_idamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_icamin
    use f77_blas, only: icamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_icamin

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = icamin(N, array, 1)
    res(2) = mfi_icamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = icamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_icamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_izamin
    use f77_blas, only: izamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_izamin

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = izamin(N, array, 1)
    res(2) = mfi_izamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = izamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_izamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_isamax
    use f77_blas, only: isamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_isamax

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = isamax(N, array, 1)
    res(2) = mfi_isamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = isamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_isamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_idamax
    use f77_blas, only: idamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_idamax

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    real(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = idamax(N, array, 1)
    res(2) = mfi_idamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(array)
    res(1) = idamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_idamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_icamax
    use f77_blas, only: icamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_icamax

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = icamax(N, array, 1)
    res(2) = mfi_icamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = icamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_icamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_izamax
    use f77_blas, only: izamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_izamax

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(wp) :: rnd(N)
    complex(wp) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = izamax(N, array, 1)
    res(2) = mfi_izamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for sequential array")

    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = izamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_izamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine

    pure subroutine assert(test, msg)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        if (.not. test) then
            error stop msg
        end if
    end subroutine
end program
