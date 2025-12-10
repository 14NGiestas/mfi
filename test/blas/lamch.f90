



program test_lamch
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

