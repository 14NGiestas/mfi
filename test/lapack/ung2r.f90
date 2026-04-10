
program test_ung2r
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cung2r 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ung2r against cung2r", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zung2r 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ung2r against zung2r", t2-t1
end block
contains

subroutine test_cung2r
    use f77_lapack, only: cung2r, f77_ung2r
    use mfi_lapack, only: mfi_ung2r, mfi_cung2r

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2
    complex(REAL32) :: A(M,N), A_copy(M,N), tau(min(M,N))
    complex(REAL32), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A in upper trapezoidal form for QR (with Householder vectors below)
    A = 0.0_wp
    A(1,1) = cmplx(2.0_wp, 0.5_wp); A(1,2) = cmplx(1.0_wp, -0.3_wp)
    A(2,1) = cmplx(1.0_wp, -0.2_wp); A(2,2) = cmplx(2.0_wp, 0.1_wp)
    A(3,1) = cmplx(1.0_wp, 0.4_wp); A(3,2) = cmplx(1.0_wp, -0.1_wp)
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call cung2r(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_cung2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cung2r")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_ung2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ung2r")

end subroutine
subroutine test_zung2r
    use f77_lapack, only: zung2r, f77_ung2r
    use mfi_lapack, only: mfi_ung2r, mfi_zung2r

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2
    complex(REAL64) :: A(M,N), A_copy(M,N), tau(min(M,N))
    complex(REAL64), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A in upper trapezoidal form for QR (with Householder vectors below)
    A = 0.0_wp
    A(1,1) = cmplx(2.0_wp, 0.5_wp); A(1,2) = cmplx(1.0_wp, -0.3_wp)
    A(2,1) = cmplx(1.0_wp, -0.2_wp); A(2,2) = cmplx(2.0_wp, 0.1_wp)
    A(3,1) = cmplx(1.0_wp, 0.4_wp); A(3,2) = cmplx(1.0_wp, -0.1_wp)
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call zung2r(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_zung2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zung2r")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    allocate(work(M))
    call mfi_ung2r(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ung2r")

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