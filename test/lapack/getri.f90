
program test_getri
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sgetri 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_getri against sgetri", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dgetri 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_getri against dgetri", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cgetri 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_getri against cgetri", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zgetri 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_getri against zgetri", t2-t1
end block
contains

subroutine test_sgetri
    use f77_lapack, only: sgetri, f77_getri, f77_getrf, sgetrf
    use mfi_lapack, only: mfi_getri, mfi_sgetri

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call sgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call sgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call sgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_sgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_sgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getri")

end subroutine
subroutine test_dgetri
    use f77_lapack, only: dgetri, f77_getri, f77_getrf, dgetrf
    use mfi_lapack, only: mfi_getri, mfi_dgetri

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call dgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call dgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call dgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_dgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_dgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getri")

end subroutine
subroutine test_cgetri
    use f77_lapack, only: cgetri, f77_getri, f77_getrf, cgetrf
    use mfi_lapack, only: mfi_getri, mfi_cgetri

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call cgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call cgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call cgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_cgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_cgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getri")

end subroutine
subroutine test_zgetri
    use f77_lapack, only: zgetri, f77_getri, f77_getrf, zgetrf
    use mfi_lapack, only: mfi_getri, mfi_zgetri

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call zgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call zgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call zgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_zgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_zgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getri")

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

    subroutine report_test_result(test_name, success)
        character(*), intent(in) :: test_name
        logical, intent(in) :: success

        if (success) then
            write(*, '(A, ": ", A)') trim(test_name), 'PASSED'
        else
            write(*, '(A, ": ", A)') trim(test_name), 'FAILED'
        end if
    end subroutine

end program
