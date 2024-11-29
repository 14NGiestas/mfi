program main
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_saxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against saxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_daxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against daxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_caxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against caxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zaxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against zaxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against sasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against dasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_scasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against scasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dzasum 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against dzasum", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_snrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against snrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dnrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against dnrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_scnrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against scnrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dznrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against dznrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_srot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against srot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_drot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against drot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_crot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against crot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against zrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_csrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against csrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against zdrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_sscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against sscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against dscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against cscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against zscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_csscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against csscal", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdscal 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing against zdscal", t2-t1
end block

contains
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

    pure subroutine assert(test, msg)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        if (.not. test) then
            error stop msg
        end if
    end subroutine
end program
