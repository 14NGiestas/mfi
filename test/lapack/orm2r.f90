
program test_orm2r
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sorm2r 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orm2r against sorm2r", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dorm2r 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orm2r against dorm2r", t2-t1
end block
contains

subroutine test_sorm2r
    use f77_lapack, only: sorm2r, f77_orm2r
    use mfi_lapack, only: mfi_orm2r, mfi_sorm2r

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2, K = min(M,N)
    real(REAL32) :: A(M,N), A_copy(M,N), C(M,N), C_orig(M,N), C_rf(M,N), tau(K)
    real(REAL32), allocatable :: work(:)
    integer :: info, info_rf, info_mfi
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A (pre-factorized form for QR)
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp
    A(2,1) = 1.0_wp; A(2,2) = 2.0_wp
    A(3,1) = 1.0_wp; A(3,2) = 1.0_wp
    
    ! Initialize tau
    tau = 0.0_wp
    tau(1) = 1.0_wp
    if (size(tau) > 1) tau(2) = 0.5_wp

    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, 1.0_wp], [M,N], order=[2,1])
    
    ! Store original for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface - need to allocate workspace
    allocate(work(M*N)) ! Workspace for multiplication routines
    call sorm2r(side, trans, M, N, K, A, M, tau, C_rf, M, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy since f77 might modify A
    allocate(work(M*N))
    call mfi_sorm2r(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sorm2r")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy since f77 might modify A
    allocate(work(M*N))
    call mfi_orm2r(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orm2r")

end subroutine
subroutine test_dorm2r
    use f77_lapack, only: dorm2r, f77_orm2r
    use mfi_lapack, only: mfi_orm2r, mfi_dorm2r

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2, K = min(M,N)
    real(REAL64) :: A(M,N), A_copy(M,N), C(M,N), C_orig(M,N), C_rf(M,N), tau(K)
    real(REAL64), allocatable :: work(:)
    integer :: info, info_rf, info_mfi
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A (pre-factorized form for QR)
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp
    A(2,1) = 1.0_wp; A(2,2) = 2.0_wp
    A(3,1) = 1.0_wp; A(3,2) = 1.0_wp
    
    ! Initialize tau
    tau = 0.0_wp
    tau(1) = 1.0_wp
    if (size(tau) > 1) tau(2) = 0.5_wp

    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, 1.0_wp], [M,N], order=[2,1])
    
    ! Store original for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface - need to allocate workspace
    allocate(work(M*N)) ! Workspace for multiplication routines
    call dorm2r(side, trans, M, N, K, A, M, tau, C_rf, M, work, info_rf)
    deallocate(work)

    ! Test mfi interface (short form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy since f77 might modify A
    allocate(work(M*N))
    call mfi_dorm2r(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dorm2r")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Use copy since f77 might modify A
    allocate(work(M*N))
    call mfi_orm2r(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    deallocate(work)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orm2r")

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