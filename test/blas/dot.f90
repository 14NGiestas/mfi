



program test_dot
use iso_fortran_env
implicit none
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
contains
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
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "different results")

    res = mfi_sdot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "different results")

    res = mfi_dot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "different results")

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
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "different results")

    res = mfi_ddot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "different results")

    res = mfi_dot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "different results")

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

