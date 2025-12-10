



program test_dotu
use iso_fortran_env
implicit none
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
contains
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

subroutine report_test_result(test_name, success)
    character(*), intent(in) :: test_name
    logical, intent(in) :: success

    if (success) then
        write(*, '(A, ": ", A)') trim(test_name), 'PASSED'
    else
        write(*, '(A, ": ", A)') trim(test_name), 'FAILED'
    end if
end subroutine

end program

