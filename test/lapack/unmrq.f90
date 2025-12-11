
program test_unmrq
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cunmrq 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_unmrq against cunmrq", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zunmrq 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_unmrq against zunmrq", t2-t1
end block
contains

subroutine test_cunmrq
    use f77_lapack, only: cunmrq, f77_unmrq
    use mfi_lapack, only: mfi_unmrq, mfi_cunmrq, mfi_geqrf, mfi_gerqf

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2
    complex(REAL32) :: A(M,M), A_copy(M,M), C(M,N), C_orig(M,N), C_rf(M,N)
    complex(REAL32) :: tau(M)
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A (will be factorized)
    ! For RQ/QR routines, the factorization matrix A needs appropriate dimensions for multiplication
    ! When using ORMRQ/ORMQR, A should be properly sized for the side operation:
    ! If side='L': A should be KxM where K=min(M,N) and we compute Q*C
    ! If side='R': A should be KxN where K=min(M,N) and we compute C*Q

    ! For the test, we'll start with an appropriately sized square matrix
    ! that can be used for both QR and RQ factorization
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 2.0_wp], [M,M])

    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, &
                 1.0_wp, 1.0_wp], [M,N])

    ! For RQ routines, factorization matrix dimensions need special handling
    ! According to documentation: if we want to use ORMRQ to multiply matrix C (M×N)
    ! A must have dimensions (LDA,M) if SIDE='L' or (LDA,N) if SIDE='R'
    ! where LDA >= max(1,K) and K follows the rule:
    ! If SIDE = 'L', M >= K >= 0; if SIDE = 'R', N >= K >= 0
    A_copy = A
    ! Compute RQ factorization first for RQ routines
    call mfi_gerqf(A_copy, tau, info=info)
    if (info /= 0) return

    ! Store original matrices for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface for cunmrq
    allocate(work(1))
    lwork = -1
    ! For RQ routines:
    ! k=min(M,N) is number of reflectors, but A matrix dimensions depend on side
    ! LDA = M if A is M×M (square matrix after factorization)
    ! For square A after RQ factorization: A is M×M, LDA=M, k=min(M,N)
    ! When side='L': A should conceptually be K×M, when side='R': A should be K×N
    ! But the factorized matrix A after GERQF is still M×M
    call cunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_copy = A  ! Reset A_copy
        call mfi_gerqf(A_copy, tau, info=info)
        if (info /= 0) then
            deallocate(work)
            return
        end if

        C_rf = C_orig  ! Reset C_rf
        call cunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form)
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_cunmrq(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cunmrq")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_unmrq(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_unmrq")

end subroutine
subroutine test_zunmrq
    use f77_lapack, only: zunmrq, f77_unmrq
    use mfi_lapack, only: mfi_unmrq, mfi_zunmrq, mfi_geqrf, mfi_gerqf

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2
    complex(REAL64) :: A(M,M), A_copy(M,M), C(M,N), C_orig(M,N), C_rf(M,N)
    complex(REAL64) :: tau(M)
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork
    character :: side = 'L', trans = 'N'

    ! Create a test matrix A (will be factorized)
    ! For RQ/QR routines, the factorization matrix A needs appropriate dimensions for multiplication
    ! When using ORMRQ/ORMQR, A should be properly sized for the side operation:
    ! If side='L': A should be KxM where K=min(M,N) and we compute Q*C
    ! If side='R': A should be KxN where K=min(M,N) and we compute C*Q

    ! For the test, we'll start with an appropriately sized square matrix
    ! that can be used for both QR and RQ factorization
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 2.0_wp], [M,M])

    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, &
                 1.0_wp, 1.0_wp], [M,N])

    ! For RQ routines, factorization matrix dimensions need special handling
    ! According to documentation: if we want to use ORMRQ to multiply matrix C (M×N)
    ! A must have dimensions (LDA,M) if SIDE='L' or (LDA,N) if SIDE='R'
    ! where LDA >= max(1,K) and K follows the rule:
    ! If SIDE = 'L', M >= K >= 0; if SIDE = 'R', N >= K >= 0
    A_copy = A
    ! Compute RQ factorization first for RQ routines
    call mfi_gerqf(A_copy, tau, info=info)
    if (info /= 0) return

    ! Store original matrices for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface for zunmrq
    allocate(work(1))
    lwork = -1
    ! For RQ routines:
    ! k=min(M,N) is number of reflectors, but A matrix dimensions depend on side
    ! LDA = M if A is M×M (square matrix after factorization)
    ! For square A after RQ factorization: A is M×M, LDA=M, k=min(M,N)
    ! When side='L': A should conceptually be K×M, when side='R': A should be K×N
    ! But the factorized matrix A after GERQF is still M×M
    call zunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_copy = A  ! Reset A_copy
        call mfi_gerqf(A_copy, tau, info=info)
        if (info /= 0) then
            deallocate(work)
            return
        end if

        C_rf = C_orig  ! Reset C_rf
        call zunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form)
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_zunmrq(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zunmrq")

    ! Test mfi interface (full form)
    C = C_orig  ! Reset C
    A_copy = A  ! Reset A_copy
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_unmrq(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_unmrq")

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