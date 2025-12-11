
program test_ormqr
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sormqr 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ormqr against sormqr", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dormqr 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ormqr against dormqr", t2-t1
end block
contains

subroutine test_sormqr
    use f77_lapack, only: sormqr, f77_ormqr
    use mfi_lapack, only: mfi_ormqr, mfi_sormqr, mfi_geqrf

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2
    real(REAL32) :: A(M,M), A_copy(M,M), C(M,N), C_orig(M,N), C_rf(M,N)
    real(REAL32) :: tau(M)
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A (will be factorized)
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 2.0_wp], [M,M])
    
    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, &
                 1.0_wp, 1.0_wp], [M,N])

    ! Compute QR factorization first
    A_copy = A
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return

    ! Store original matrices for comparison
    C_orig = C
    C_rf = C
    
    ! Test f77 interface for ormqr
    allocate(work(1))
    lwork = -1
    call sormqr(side, trans, M, N, M, A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_copy = A  ! Reset A_copy
        call mfi_geqrf(A_copy, tau, info=info)
        if (info /= 0) then
            deallocate(work)
            return
        end if
        
        C_rf = C_orig  ! Reset C_rf
        call sormqr(side, trans, M, N, M, A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form) - this should apply Q*C operation
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_sormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sormqr")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_ormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ormqr")

end subroutine
subroutine test_dormqr
    use f77_lapack, only: dormqr, f77_ormqr
    use mfi_lapack, only: mfi_ormqr, mfi_dormqr, mfi_geqrf

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2
    real(REAL64) :: A(M,M), A_copy(M,M), C(M,N), C_orig(M,N), C_rf(M,N)
    real(REAL64) :: tau(M)
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A (will be factorized)
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 2.0_wp], [M,M])
    
    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, &
                 1.0_wp, 1.0_wp], [M,N])

    ! Compute QR factorization first
    A_copy = A
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return

    ! Store original matrices for comparison
    C_orig = C
    C_rf = C
    
    ! Test f77 interface for ormqr
    allocate(work(1))
    lwork = -1
    call dormqr(side, trans, M, N, M, A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_copy = A  ! Reset A_copy
        call mfi_geqrf(A_copy, tau, info=info)
        if (info /= 0) then
            deallocate(work)
            return
        end if
        
        C_rf = C_orig  ! Reset C_rf
        call dormqr(side, trans, M, N, M, A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form) - this should apply Q*C operation
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_dormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dormqr")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_ormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ormqr")

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