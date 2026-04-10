
program test_orgr2
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sorgr2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orgr2 against sorgr2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dorgr2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orgr2 against dorgr2", t2-t1
end block
contains

subroutine test_sorgr2
    use f77_lapack, only: sorgr2, f77_orgr2
    use mfi_lapack, only: mfi_orgr2, mfi_sorgr2

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 2, N = 3  ! For RQ, typically L <= N, so M <= N
    real(REAL32) :: A(M,N), A_copy(M,N), tau(min(M,N))
    real(REAL32), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A in lower trapezoidal form for RQ (with Householder vectors above)
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp; A(1,3) = 1.0_wp
    A(2,1) = 1.0_wp; A(2,2) = 2.0_wp; A(2,3) = 1.0_wp
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call sorgr2(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_sorgr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sorgr2")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_orgr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgr2")

end subroutine
subroutine test_dorgr2
    use f77_lapack, only: dorgr2, f77_orgr2
    use mfi_lapack, only: mfi_orgr2, mfi_dorgr2

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 2, N = 3  ! For RQ, typically L <= N, so M <= N
    real(REAL64) :: A(M,N), A_copy(M,N), tau(min(M,N))
    real(REAL64), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A in lower trapezoidal form for RQ (with Householder vectors above)
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp; A(1,3) = 1.0_wp
    A(2,1) = 1.0_wp; A(2,2) = 2.0_wp; A(2,3) = 1.0_wp
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call dorgr2(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_dorgr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dorgr2")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_orgr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgr2")

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