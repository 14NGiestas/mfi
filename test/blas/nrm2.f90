

program test_nrm2
use iso_fortran_env
implicit none
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
contains
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

