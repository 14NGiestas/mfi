
program test_ungr2
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cungr2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ungr2 against cungr2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zungr2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ungr2 against zungr2", t2-t1
end block
contains

subroutine test_cungr2
    use f77_lapack, only: cungr2, f77_ungr2
    use mfi_lapack, only: mfi_ungr2, mfi_cungr2

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 2, N = 3
    complex(REAL32) :: A(M,N), A_copy(M,N), tau(min(M,N))
    complex(REAL32), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A 
    A = 0.0_wp
    A(1,1) = cmplx(2.0_wp, 0.5_wp); A(1,2) = cmplx(1.0_wp, -0.3_wp); A(1,3) = cmplx(1.0_wp, 0.2_wp)
    A(2,1) = cmplx(1.0_wp, -0.2_wp); A(2,2) = cmplx(2.0_wp, 0.1_wp); A(2,3) = cmplx(1.0_wp, -0.1_wp)
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call cungr2(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A
    allocate(work(M))
    call mfi_cungr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cungr2")

    ! Test mfi interface (full form)
    A_copy = A
    allocate(work(M))
    call mfi_ungr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungr2")

end subroutine
subroutine test_zungr2
    use f77_lapack, only: zungr2, f77_ungr2
    use mfi_lapack, only: mfi_ungr2, mfi_zungr2

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 2, N = 3
    complex(REAL64) :: A(M,N), A_copy(M,N), tau(min(M,N))
    complex(REAL64), allocatable :: work(:)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix A 
    A = 0.0_wp
    A(1,1) = cmplx(2.0_wp, 0.5_wp); A(1,2) = cmplx(1.0_wp, -0.3_wp); A(1,3) = cmplx(1.0_wp, 0.2_wp)
    A(2,1) = cmplx(1.0_wp, -0.2_wp); A(2,2) = cmplx(2.0_wp, 0.1_wp); A(2,3) = cmplx(1.0_wp, -0.1_wp)
    
    A_copy = A

    ! Test f77 interface - need to allocate workspace
    allocate(work(M)) ! Workspace size for generation routines
    call zungr2(M, N, min(M,N), A_copy, M, tau, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    A_copy = A
    allocate(work(M))
    call mfi_zungr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zungr2")

    ! Test mfi interface (full form)
    A_copy = A
    allocate(work(M))
    call mfi_ungr2(A_copy, tau, k=min(M,N), info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(A - A_copy) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungr2")

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