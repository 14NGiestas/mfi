
program test_asum
use iso_fortran_env
implicit none
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
contains
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
