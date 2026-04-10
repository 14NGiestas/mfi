
program test_unmr2
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cunmr2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_unmr2 against cunmr2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zunmr2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_unmr2 against zunmr2", t2-t1
end block
contains

subroutine test_cunmr2
    use f77_lapack, only: cunmr2, f77_unmr2
    use mfi_lapack, only: mfi_unmr2, mfi_cunmr2

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2, K = min(M,N)
    complex(REAL32) :: A(M,N), A_copy(M,N), C(M,N), C_orig(M,N), C_rf(M,N), tau(K)
    complex(REAL32), allocatable :: work(:)
    integer :: info, info_rf, info_mfi
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A in appropriate form for RQ factorization
    A = 0.0_wp
    A(1,1) = cmplx(2.0_wp, 0.5_wp); A(1,2) = cmplx(1.0_wp, -0.3_wp)
    A(2,1) = cmplx(1.0_wp, -0.2_wp); A(2,2) = cmplx(2.0_wp, 0.1_wp)
    A(3,1) = cmplx(1.0_wp, 0.4_wp); A(3,2) = cmplx(1.0_wp, -0.1_wp)
    
    ! Initialize tau for RQ
    tau = cmplx(0.0_wp, 0.0_wp)
    tau(1) = cmplx(1.0_wp, 0.5_wp)
    if (size(tau) > 1) tau(2) = cmplx(0.5_wp, -0.2_wp)

    ! Create test matrix C for multiplication
    C = reshape([cmplx(2.0_wp, 0.1_wp), cmplx(1.0_wp, -0.1_wp), cmplx(1.0_wp, 0.2_wp), &
                 cmplx(1.0_wp, -0.1_wp), cmplx(2.0_wp, 0.3_wp), cmplx(1.0_wp, -0.2_wp)], [M,N], order=[2,1])
    
    ! Store original for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface - need to allocate workspace
    allocate(work(M*N)) ! Workspace for multiplication routines
    call cunmr2(side, trans, M, N, K, A, M, tau, C_rf, M, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy
    allocate(work(M*N))
    call mfi_cunmr2(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cunmr2")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy
    allocate(work(M*N))
    call mfi_unmr2(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_unmr2")

end subroutine
subroutine test_zunmr2
    use f77_lapack, only: zunmr2, f77_unmr2
    use mfi_lapack, only: mfi_unmr2, mfi_zunmr2

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2, K = min(M,N)
    complex(REAL64) :: A(M,N), A_copy(M,N), C(M,N), C_orig(M,N), C_rf(M,N), tau(K)
    complex(REAL64), allocatable :: work(:)
    integer :: info, info_rf, info_mfi
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A in appropriate form for RQ factorization
    A = 0.0_wp
    A(1,1) = cmplx(2.0_wp, 0.5_wp); A(1,2) = cmplx(1.0_wp, -0.3_wp)
    A(2,1) = cmplx(1.0_wp, -0.2_wp); A(2,2) = cmplx(2.0_wp, 0.1_wp)
    A(3,1) = cmplx(1.0_wp, 0.4_wp); A(3,2) = cmplx(1.0_wp, -0.1_wp)
    
    ! Initialize tau for RQ
    tau = cmplx(0.0_wp, 0.0_wp)
    tau(1) = cmplx(1.0_wp, 0.5_wp)
    if (size(tau) > 1) tau(2) = cmplx(0.5_wp, -0.2_wp)

    ! Create test matrix C for multiplication
    C = reshape([cmplx(2.0_wp, 0.1_wp), cmplx(1.0_wp, -0.1_wp), cmplx(1.0_wp, 0.2_wp), &
                 cmplx(1.0_wp, -0.1_wp), cmplx(2.0_wp, 0.3_wp), cmplx(1.0_wp, -0.2_wp)], [M,N], order=[2,1])
    
    ! Store original for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface - need to allocate workspace
    allocate(work(M*N)) ! Workspace for multiplication routines
    call zunmr2(side, trans, M, N, K, A, M, tau, C_rf, M, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy
    allocate(work(M*N))
    call mfi_zunmr2(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zunmr2")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy
    allocate(work(M*N))
    call mfi_unmr2(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_unmr2")

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