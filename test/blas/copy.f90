



program test_copy
use iso_fortran_env
implicit none
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
contains
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
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_scopy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

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
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_dcopy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

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
    x = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im, kind=REAL32)
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
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_ccopy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

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
    x = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im, kind=REAL64)
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
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_zcopy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = x
    y_in = y
    call mfi_copy(x_in, y_in)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

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

