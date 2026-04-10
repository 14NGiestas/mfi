



program test_scal
use iso_fortran_env
implicit none
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
contains
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
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_sscal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

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
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_dscal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

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
    X = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL32)
end block

    ! The test is always against the original
    x_in = x
    call cscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_cscal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

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
    X = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL64)
end block

    ! The test is always against the original
    x_in = x
    call zscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_zscal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

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
    X = cmplx(re,im, kind=REAL32)
end block
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call csscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_csscal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

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
    X = cmplx(re,im, kind=REAL64)
end block
    call random_number(alpha)

    ! The test is always against the original
    x_in = x
    call zdscal(N, alpha, x_in, 1)
    x_rf = x_in

    x_in = x
    call f77_scal(N, alpha, x_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_zdscal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

    x_in = x
    call mfi_scal(alpha, x_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))), "different results")

end subroutine

subroutine assert(test, msg, info)
    logical, intent(in) :: test
    character(*), intent(in) :: msg
    integer, intent(in), optional :: info
    character(1024) :: buffer

    if (.not. test) then
        if (present(info)) then
            write(buffer, *) 'Error ', info, ': ', msg
        else
            write(buffer, *) 'Error: ', msg
        end if
        error stop trim(buffer)
    end if
end subroutine

end program

