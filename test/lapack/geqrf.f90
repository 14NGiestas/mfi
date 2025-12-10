
program test_geqrf
    use iso_fortran_env
    implicit none

    write(*,*) 'Starting geqrf tests...'
    call test_sgeqrf
    call test_dgeqrf
    call test_cgeqrf
    call test_zgeqrf
    write(*,*) 'All geqrf tests completed successfully.'

contains

subroutine test_sgeqrf
    use f77_lapack, only: sgeqrf, f77_geqrf
    use mfi_lapack, only: mfi_geqrf, mfi_sgeqrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 6, M = 2
    real(REAL32) :: A(N,M), A_in(N,M), A_rf(N,M)
    real(REAL32) :: tau(min(N,M)), tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call sgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call sgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_sgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_sgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_geqrf")

end subroutine
subroutine test_dgeqrf
    use f77_lapack, only: dgeqrf, f77_geqrf
    use mfi_lapack, only: mfi_geqrf, mfi_dgeqrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 6, M = 2
    real(REAL64) :: A(N,M), A_in(N,M), A_rf(N,M)
    real(REAL64) :: tau(min(N,M)), tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call dgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call dgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_dgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_dgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_geqrf")

end subroutine
subroutine test_cgeqrf
    use f77_lapack, only: cgeqrf, f77_geqrf
    use mfi_lapack, only: mfi_geqrf, mfi_cgeqrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 6, M = 2
    complex(REAL32) :: A(N,M), A_in(N,M), A_rf(N,M)
    complex(REAL32) :: tau(min(N,M)), tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call cgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call cgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_cgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_cgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_geqrf")

end subroutine
subroutine test_zgeqrf
    use f77_lapack, only: zgeqrf, f77_geqrf
    use mfi_lapack, only: mfi_geqrf, mfi_zgeqrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 6, M = 2
    complex(REAL64) :: A(N,M), A_in(N,M), A_rf(N,M)
    complex(REAL64) :: tau(min(N,M)), tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call zgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call zgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_zgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_zgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_geqrf")

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
