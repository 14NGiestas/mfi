program main
use iso_fortran_env
implicit none
print*, 'calling sasum'
call test_sasum
print*, 'calling dasum'
call test_dasum
print*, 'calling scasum'
call test_scasum
print*, 'calling dzasum'
call test_dzasum

print*, 'calling snrm2'
call test_snrm2
print*, 'calling dnrm2'
call test_dnrm2
print*, 'calling scnrm2'
call test_scnrm2
print*, 'calling dznrm2'
call test_dznrm2
contains
subroutine test_sasum
    use f77_blas, only: sasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_sasum

    integer, parameter :: wp = REAL32 
    integer, parameter :: N = 20
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
    res(2) = mfi_sasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dasum
    use f77_blas, only: dasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_dasum

    integer, parameter :: wp = REAL64 
    integer, parameter :: N = 20
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
    res(2) = mfi_dasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_scasum
    use f77_blas, only: scasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_scasum

    integer, parameter :: wp = REAL32 
    integer, parameter :: N = 20
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

    call random_number(array%re)
    call random_number(array%im)
    res(1) = scasum(N, array, 1)
    res(2) = mfi_scasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dzasum
    use f77_blas, only: dzasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_dzasum

    integer, parameter :: wp = REAL64 
    integer, parameter :: N = 20
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

    call random_number(array%re)
    call random_number(array%im)
    res(1) = dzasum(N, array, 1)
    res(2) = mfi_dzasum(array)
    res(3) = f77_asum(N, array, 1)
    res(4) = mfi_asum(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine

subroutine test_snrm2
    use f77_blas, only: snrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_snrm2

    integer, parameter :: wp = REAL32 
    integer, parameter :: N = 20
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
    res(2) = mfi_snrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dnrm2
    use f77_blas, only: dnrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_dnrm2

    integer, parameter :: wp = REAL64 
    integer, parameter :: N = 20
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
    res(2) = mfi_dnrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_scnrm2
    use f77_blas, only: scnrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_scnrm2

    integer, parameter :: wp = REAL32 
    integer, parameter :: N = 20
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

    call random_number(array%re)
    call random_number(array%im)
    res(1) = scnrm2(N, array, 1)
    res(2) = mfi_scnrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
    call assert(all(res == res(1)), "different results for random array")

end subroutine
subroutine test_dznrm2
    use f77_blas, only: dznrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_dznrm2

    integer, parameter :: wp = REAL64 
    integer, parameter :: N = 20
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

    call random_number(array%re)
    call random_number(array%im)
    res(1) = dznrm2(N, array, 1)
    res(2) = mfi_dznrm2(array)
    res(3) = f77_nrm2(N, array, 1)
    res(4) = mfi_nrm2(array)
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
