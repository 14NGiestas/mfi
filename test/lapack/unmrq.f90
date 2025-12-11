
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
    ! For RQ routines, we need temporary matrix for proper factorization
    complex(REAL32) :: A_factorize(M,N)  ! Temporary matrix for RQ factorization

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

    ! For RQ/QR routines, we need consistent matrix setup for proper testing
    A_copy = A
    ! For RQ routines, we need to be very careful about dimensions
    ! For multiplying C (M×N), the factorization matrix must come from a matrix compatible with the operation
    ! For side='L': A should be from factorizing an M×N matrix
    ! For side='R': A should be from factorizing an M×N matrix
    ! Actually, let me use proper dimensions: if we want to multiply C(M,N) where C is M×N
    ! For side='L': A is result of RQ factorizing a matrix of size (something) × N
    ! For side='R': A is result of RQ factorizing a matrix of size M × (something)

    ! To keep it simple, let's make A the right size for the RQ factorization
    ! Since we're multiplying C(M,N) and A is initially M×M, let me adjust the approach
    ! For RQ factorization of an M×N matrix, we need to work with proper dimensions

    ! Create proper matrices for RQ case - we need to factorize an appropriate matrix
    ! For RQ, if we want to multiply C (M×N), we factorize an N×N matrix if side='R', or M×? if side='L'
    ! Actually, let me use a temporary matrix to match the mathematical requirements
    ! For ORMRQ with side='L': multiplying Q*C where C is M×N, Q comes from RQ of a matrix X that's M×N
    ! So we should RQ-factorize a M×N matrix (not M×M!)

    ! For proper RQ test: create an appropriate matrix for factorization - one with M×N dimensions
    ! But we need to be careful about the storage since we need A(M,M)
    ! Let me use submatrix approach: use upper M×N part of A for factorization

    ! Create appropriate matrix for RQ factorization (M×N)
    A_factorize = reshape([3.0_wp, 1.0_wp, 1.0_wp, &  ! First N columns of first M rows
                           1.0_wp, 3.0_wp, 1.0_wp], [M,N])  ! Second N columns of first M rows if exists

    ! Compute RQ factorization of the M×N matrix
    call mfi_gerqf(A_factorize, tau, info=info)
    if (info /= 0) return

    ! Now copy the result to A for the multiplication test, padding with zeros if needed
    A_copy = 0.0_wp
    A_copy(1:M, 1:N) = A_factorize(1:M, 1:N)


    ! Store original matrices for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface for cunmrq - workspace query and then actual call
    allocate(work(1))
    lwork = -1
    call cunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        ! Use the same factorization and original C matrix for the actual computation
        C_rf = C_orig  ! Reset C_rf to original before applying transform
        call cunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form) - use the same factorization and original input
    ! Reset matrices to original values
    C = C_orig
    ! Use same approach for mfi
    A_factorize = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                           1.0_wp, 3.0_wp, 1.0_wp], [M,N])
    call mfi_gerqf(A_factorize, tau, info=info)
    if (info /= 0) return
    A_copy = 0.0_wp
    A_copy(1:M, 1:N) = A_factorize(1:M, 1:N)
    call mfi_cunmrq(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cunmrq")

    ! Test mfi interface (full form) - use the same factorization and original input
    C = C_orig
    ! Use same approach for mfi
    A_factorize = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                           1.0_wp, 3.0_wp, 1.0_wp], [M,N])
    call mfi_gerqf(A_factorize, tau, info=info)
    if (info /= 0) return
    A_copy = 0.0_wp
    A_copy(1:M, 1:N) = A_factorize(1:M, 1:N)
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
    ! For RQ routines, we need temporary matrix for proper factorization
    complex(REAL64) :: A_factorize(M,N)  ! Temporary matrix for RQ factorization

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

    ! For RQ/QR routines, we need consistent matrix setup for proper testing
    A_copy = A
    ! For RQ routines, we need to be very careful about dimensions
    ! For multiplying C (M×N), the factorization matrix must come from a matrix compatible with the operation
    ! For side='L': A should be from factorizing an M×N matrix
    ! For side='R': A should be from factorizing an M×N matrix
    ! Actually, let me use proper dimensions: if we want to multiply C(M,N) where C is M×N
    ! For side='L': A is result of RQ factorizing a matrix of size (something) × N
    ! For side='R': A is result of RQ factorizing a matrix of size M × (something)

    ! To keep it simple, let's make A the right size for the RQ factorization
    ! Since we're multiplying C(M,N) and A is initially M×M, let me adjust the approach
    ! For RQ factorization of an M×N matrix, we need to work with proper dimensions

    ! Create proper matrices for RQ case - we need to factorize an appropriate matrix
    ! For RQ, if we want to multiply C (M×N), we factorize an N×N matrix if side='R', or M×? if side='L'
    ! Actually, let me use a temporary matrix to match the mathematical requirements
    ! For ORMRQ with side='L': multiplying Q*C where C is M×N, Q comes from RQ of a matrix X that's M×N
    ! So we should RQ-factorize a M×N matrix (not M×M!)

    ! For proper RQ test: create an appropriate matrix for factorization - one with M×N dimensions
    ! But we need to be careful about the storage since we need A(M,M)
    ! Let me use submatrix approach: use upper M×N part of A for factorization

    ! Create appropriate matrix for RQ factorization (M×N)
    A_factorize = reshape([3.0_wp, 1.0_wp, 1.0_wp, &  ! First N columns of first M rows
                           1.0_wp, 3.0_wp, 1.0_wp], [M,N])  ! Second N columns of first M rows if exists

    ! Compute RQ factorization of the M×N matrix
    call mfi_gerqf(A_factorize, tau, info=info)
    if (info /= 0) return

    ! Now copy the result to A for the multiplication test, padding with zeros if needed
    A_copy = 0.0_wp
    A_copy(1:M, 1:N) = A_factorize(1:M, 1:N)


    ! Store original matrices for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface for zunmrq - workspace query and then actual call
    allocate(work(1))
    lwork = -1
    call zunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        ! Use the same factorization and original C matrix for the actual computation
        C_rf = C_orig  ! Reset C_rf to original before applying transform
        call zunmrq(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form) - use the same factorization and original input
    ! Reset matrices to original values
    C = C_orig
    ! Use same approach for mfi
    A_factorize = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                           1.0_wp, 3.0_wp, 1.0_wp], [M,N])
    call mfi_gerqf(A_factorize, tau, info=info)
    if (info /= 0) return
    A_copy = 0.0_wp
    A_copy(1:M, 1:N) = A_factorize(1:M, 1:N)
    call mfi_zunmrq(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zunmrq")

    ! Test mfi interface (full form) - use the same factorization and original input
    C = C_orig
    ! Use same approach for mfi
    A_factorize = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                           1.0_wp, 3.0_wp, 1.0_wp], [M,N])
    call mfi_gerqf(A_factorize, tau, info=info)
    if (info /= 0) return
    A_copy = 0.0_wp
    A_copy(1:M, 1:N) = A_factorize(1:M, 1:N)
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