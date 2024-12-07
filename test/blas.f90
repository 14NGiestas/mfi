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
 call test_sdot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dot against sdot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ddot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dot against ddot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cdotc 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dotc against cdotc", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdotc 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dotc against zdotc", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cdotu 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dotu against cdotu", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdotu 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dotu against zdotu", t2-t1
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
 call test_srotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against srotg", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_drotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against drotg", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_crotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against crotg", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zrotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against zrotg", t2-t1
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
    real(REAL32) :: a, b
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
    real(REAL64) :: a, b
    integer :: i

    do i=1,size(options)
        a = dlamch(options(i))
        b = mfi_lamch(options(i),1.0_wp)
        call assert(a == b, "different results for option "//options(i))
    end do

end subroutine
subroutine test_sdot
    use f77_blas, only: sdot, f77_dot
    use mfi_blas, only: mfi_dot, mfi_sdot

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(REAL32) :: res, ref

    real(REAL32) :: x(N), y(N)

    call random_number(X)
    call random_number(Y)

    ! The test is always against the original
    ref = sdot(N, x, 1, y, 1)

    res = f77_dot(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_sdot(x, y)
    call assert(ref == res, "different results")

    res = mfi_dot(x, y)
    call assert(ref == res, "different results")

end subroutine
subroutine test_ddot
    use f77_blas, only: ddot, f77_dot
    use mfi_blas, only: mfi_dot, mfi_ddot

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(REAL64) :: res, ref

    real(REAL64) :: x(N), y(N)

    call random_number(X)
    call random_number(Y)

    ! The test is always against the original
    ref = ddot(N, x, 1, y, 1)

    res = f77_dot(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_ddot(x, y)
    call assert(ref == res, "different results")

    res = mfi_dot(x, y)
    call assert(ref == res, "different results")

end subroutine
subroutine test_cdotc
    use f77_blas, only: cdotc, f77_dotc
    use mfi_blas, only: mfi_dotc, mfi_cdotc

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    complex(REAL32) :: res, ref

    complex(REAL32) :: x(N), y(N)

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

    ! The test is always against the original
    ref = cdotc(N, x, 1, y, 1)

    res = f77_dotc(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_cdotc(x, y)
    call assert(ref == res, "different results")

    res = mfi_dotc(x, y)
    call assert(ref == res, "different results")

end subroutine
subroutine test_zdotc
    use f77_blas, only: zdotc, f77_dotc
    use mfi_blas, only: mfi_dotc, mfi_zdotc

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    complex(REAL64) :: res, ref

    complex(REAL64) :: x(N), y(N)

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

    ! The test is always against the original
    ref = zdotc(N, x, 1, y, 1)

    res = f77_dotc(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_zdotc(x, y)
    call assert(ref == res, "different results")

    res = mfi_dotc(x, y)
    call assert(ref == res, "different results")

end subroutine
subroutine test_cdotu
    use f77_blas, only: cdotu, f77_dotu
    use mfi_blas, only: mfi_dotu, mfi_cdotu

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    complex(REAL32) :: res, ref

    complex(REAL32) :: x(N), y(N)

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

    ! The test is always against the original
    ref = cdotu(N, x, 1, y, 1)

    res = f77_dotu(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_cdotu(x, y)
    call assert(ref == res, "different results")

    res = mfi_dotu(x, y)
    call assert(ref == res, "different results")

end subroutine
subroutine test_zdotu
    use f77_blas, only: zdotu, f77_dotu
    use mfi_blas, only: mfi_dotu, mfi_zdotu

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    complex(REAL64) :: res, ref

    complex(REAL64) :: x(N), y(N)

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

    ! The test is always against the original
    ref = zdotu(N, x, 1, y, 1)

    res = f77_dotu(N, x, 1, y, 1)
    call assert(ref == res, "different results")

    res = mfi_zdotu(x, y)
    call assert(ref == res, "different results")

    res = mfi_dotu(x, y)
    call assert(ref == res, "different results")

end subroutine
subroutine test_scopy
    use f77_blas, only: scopy, f77_copy
    use mfi_blas, only: mfi_copy, mfi_scopy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(REAL32) :: x(N),    y(N),    &
                    x_in(N), y_in(N), &
                    x_rf(N), y_rf(N)

    call random_number(x)
    call random_number(y)

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

    real(REAL64) :: x(N),    y(N),    &
                    x_in(N), y_in(N), &
                    x_rf(N), y_rf(N)

    call random_number(x)
    call random_number(y)

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

    complex(REAL32) :: x(N),    y(N),    &
                    x_in(N), y_in(N), &
                    x_rf(N), y_rf(N)

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    x = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im)
end block

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

    complex(REAL64) :: x(N),    y(N),    &
                    x_in(N), y_in(N), &
                    x_rf(N), y_rf(N)

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    x = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im)
end block

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

    real(REAL32) :: rnd(N)

    real(REAL32) :: x(N),    y(N),    &
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

    real(REAL64) :: rnd(N)

    real(REAL64) :: x(N),    y(N),    &
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

    real(REAL32) :: rnd(N)

    complex(REAL32) :: x(N),    y(N),    &
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

    real(REAL64) :: rnd(N)

    complex(REAL64) :: x(N),    y(N),    &
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
    real(REAL32) :: rnd_vector(N), rnd
    real(REAL32) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(REAL32) :: alpha

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
    real(REAL64) :: rnd_vector(N), rnd
    real(REAL64) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(REAL64) :: alpha

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
    real(REAL32) :: rnd_vector(N), rnd
    complex(REAL32) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    complex(REAL32) :: alpha

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
    real(REAL64) :: rnd_vector(N), rnd
    complex(REAL64) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    complex(REAL64) :: alpha

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
    real(REAL32) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    real(REAL32) :: alpha, beta
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
    real(REAL64) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    real(REAL64) :: alpha, beta
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
    complex(REAL32) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    complex(REAL32) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
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
    complex(REAL64) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    complex(REAL64) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
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
    real(REAL32) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    real(REAL32) :: alpha, beta
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
    real(REAL64) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    real(REAL64) :: alpha, beta
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
    complex(REAL32) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    complex(REAL32) :: alpha, beta
    character :: transa, transb
    integer :: i, j

block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im)
end block
block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im)
end block
block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    C = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
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
    complex(REAL64) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    complex(REAL64) :: alpha, beta
    character :: transa, transb
    integer :: i, j

block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im)
end block
block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im)
end block
block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    C = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
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
    real(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = sasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_sasum(array)
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
    real(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_dasum(array)
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
    complex(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = scasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_scasum(array)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for sequential array")

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im)
end block
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
    complex(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dzasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_dzasum(array)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for sequential array")

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im)
end block
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
    real(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = snrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_snrm2(array)
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
    real(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dnrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_dnrm2(array)
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
    complex(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = scnrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_scnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for sequential array")

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im)
end block
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
    complex(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = dznrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_dznrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for sequential array")

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im)
end block
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
    real(REAL32) :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL32) :: c
    real(REAL32) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(X)
    call random_number(Y)

    c = cos(angle)
    s = sin(angle)

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
    real(REAL64) :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL64) :: c
    real(REAL64) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

    call random_number(X)
    call random_number(Y)

    c = cos(angle)
    s = sin(angle)

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
    complex(REAL32) :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL32) :: c
    complex(REAL32) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

    c = cos(angle)
    s = sin(angle)

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
    complex(REAL64) :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL64) :: c
    complex(REAL64) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

    c = cos(angle)
    s = sin(angle)

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
    complex(REAL32) :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL32) :: c
    real(REAL32) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

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
    complex(REAL64) :: x(N),    y(N),    &
                   x_in(N), y_in(N), &
                   x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL64) :: c
    real(REAL64) :: s

    call random_number(angle)
    angle = angle * 2.0_wp * pi

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block

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
subroutine test_srotg
    use f77_blas, only: srotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 200
    real(REAL32) :: a, b, s
    real(REAL32) :: c

    real(REAL32) :: a_in, b_in, s_in
    real(REAL32) :: c_in

    real(REAL32) :: a_rf, b_rf, s_rf
    real(REAL32) :: c_rf
    integer :: i

    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(s)

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call srotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_drotg
    use f77_blas, only: drotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 200
    real(REAL64) :: a, b, s
    real(REAL64) :: c

    real(REAL64) :: a_in, b_in, s_in
    real(REAL64) :: c_in

    real(REAL64) :: a_rf, b_rf, s_rf
    real(REAL64) :: c_rf
    integer :: i

    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(s)

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call drotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_crotg
    use f77_blas, only: crotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 200
    complex(REAL32) :: a, b, s
    real(REAL32) :: c

    complex(REAL32) :: a_in, b_in, s_in
    real(REAL32) :: c_in

    complex(REAL32) :: a_rf, b_rf, s_rf
    real(REAL32) :: c_rf
    integer :: i

block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    a = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    b = cmplx(re,im)
end block
    call random_number(c)
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    s = cmplx(re,im)
end block

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call crotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_zrotg
    use f77_blas, only: zrotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 200
    complex(REAL64) :: a, b, s
    real(REAL64) :: c

    complex(REAL64) :: a_in, b_in, s_in
    real(REAL64) :: c_in

    complex(REAL64) :: a_rf, b_rf, s_rf
    real(REAL64) :: c_rf
    integer :: i

block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    a = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    b = cmplx(re,im)
end block
    call random_number(c)
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    s = cmplx(re,im)
end block

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call zrotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_sscal
    use f77_blas, only: sscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_sscal

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(REAL32) :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    real(REAL32) :: alpha

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
    call mfi_sscal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_dscal
    use f77_blas, only: dscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_dscal

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    real(REAL64) :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    real(REAL64) :: alpha

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
    call mfi_dscal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_cscal
    use f77_blas, only: cscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_cscal

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    complex(REAL32) :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    complex(REAL32) :: alpha

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block

    ! The test is always against the original
    x_in = x
    call cscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_cscal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_zscal
    use f77_blas, only: zscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_zscal

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    complex(REAL64) :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    complex(REAL64) :: alpha

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block

    ! The test is always against the original
    x_in = x
    call zscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_zscal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_csscal
    use f77_blas, only: csscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_csscal

    integer,     parameter :: wp = REAL32
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    complex(REAL32) :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    real(REAL32) :: alpha

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call csscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_csscal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine
subroutine test_zdscal
    use f77_blas, only: zdscal, f77_scal
    use mfi_blas, only: mfi_scal, mfi_zdscal

    integer,     parameter :: wp = REAL64
    real(wp),    parameter :: pi = 4.0_wp * atan(1.0_wp)
    complex(wp), parameter :: i = (0.0_wp,1.0_wp)

    integer, parameter :: N = 20
    complex(REAL64) :: x(N),    &
                   x_in(N), &
                   x_rf(N)
    real(REAL64) :: alpha

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call zdscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_zdscal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(x_in == x_rf), "different results")

end subroutine

subroutine test_isamin
    use f77_blas, only: isamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_isamin

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: rnd(N)
    real(REAL32) :: array(N)
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
    real(REAL64) :: rnd(N)
    real(REAL64) :: array(N)
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
    real(REAL32) :: rnd(N)
    complex(REAL32) :: array(N)
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
    real(REAL64) :: rnd(N)
    complex(REAL64) :: array(N)
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
    real(REAL32) :: rnd(N)
    real(REAL32) :: array(N)
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
    real(REAL64) :: rnd(N)
    real(REAL64) :: array(N)
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
    real(REAL32) :: rnd(N)
    complex(REAL32) :: array(N)
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
    real(REAL64) :: rnd(N)
    complex(REAL64) :: array(N)
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
