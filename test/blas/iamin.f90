

program test_iamin
use iso_fortran_env
implicit none
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
contains
subroutine test_isamin
    use f77_blas, only: isamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_isamin

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = isamin(N, array, 1)
    res(2) = mfi_isamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for sequential array")

    call random_number(array)
    res(1) = isamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_isamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for random array")

end subroutine
subroutine test_idamin
    use f77_blas, only: idamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_idamin

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = idamin(N, array, 1)
    res(2) = mfi_idamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for sequential array")

    call random_number(array)
    res(1) = idamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_idamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for random array")

end subroutine
subroutine test_icamin
    use f77_blas, only: icamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_icamin

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = icamin(N, array, 1)
    res(2) = mfi_icamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for sequential array")

    real(REAL32) :: rnd(N)
    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = icamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_icamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for random array")

end subroutine
subroutine test_izamin
    use f77_blas, only: izamin, f77_iamin
    use mfi_blas, only: mfi_iamin, mfi_izamin

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: array(N)
    integer :: res(4)
    integer :: i

    ! Test sequential array
    array = [(1.0_wp*i,i=1,N)]
    res(1) = izamin(N, array, 1)
    res(2) = mfi_izamin(array)
    res(3) = f77_iamin(N, array, 1)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for sequential array")

    real(REAL64) :: rnd(N)
    call random_number(rnd)
    array%re = rnd
    call random_number(rnd)
    array%im = rnd
    res(1) = izamin(N, array, 1)
    res(2) = f77_iamin(N, array, 1)
    res(3) = mfi_izamin(array)
    res(4) = mfi_iamin(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "different results for random array")

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

