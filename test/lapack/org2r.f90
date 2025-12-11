
program test_org2r
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sorg2r 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_org2r against sorg2r", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dorg2r 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_org2r against dorg2r", t2-t1
end block
contains

subroutine test_sorg2r
    use f77_lapack, only: sorg2r, f77_org2r
    use mfi_lapack, only: mfi_org2r, mfi_sorg2r

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2
    real(REAL32) :: A(M,N), A_copy(M,N), tau(min(M,N))
    real(REAL32), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A in upper trapezoidal form for QR (with Householder vectors below)
    ! For GEQR2, input should be in the form produced by the factorization
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp
    A(2,1) = 1.0_wp; A(2,2) = 2.0_wp
    A(3,1) = 1.0_wp; A(3,2) = 1.0_wp
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call sorg2r(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_sorg2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sorg2r")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_org2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_org2r")

end subroutine
subroutine test_dorg2r
    use f77_lapack, only: dorg2r, f77_org2r
    use mfi_lapack, only: mfi_org2r, mfi_dorg2r

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2
    real(REAL64) :: A(M,N), A_copy(M,N), tau(min(M,N))
    real(REAL64), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A in upper trapezoidal form for QR (with Householder vectors below)
    ! For GEQR2, input should be in the form produced by the factorization
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp
    A(2,1) = 1.0_wp; A(2,2) = 2.0_wp
    A(3,1) = 1.0_wp; A(3,2) = 1.0_wp
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call dorg2r(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_dorg2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dorg2r")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_org2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_org2r")

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