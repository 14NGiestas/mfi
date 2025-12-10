
program test_gesvd
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting gesvd tests...'
    call test_sgesvd
    call test_dgesvd
    call test_cgesvd
    call test_zgesvd
    write(*,'(A)') 'All gesvd tests completed successfully.'

contains

subroutine test_sgesvd
    use f77_lapack, only: sgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_sgesvd

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 3
    real(REAL32) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL32) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    real(REAL32), allocatable :: rwork(:)  ! Needed for complex types
    integer :: lwork

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    allocate(work(1))
    lwork = -1  ! Workspace query
        ! For real types, no rwork needed
        call sgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
            call sgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_sgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_sgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_gesvd")

end subroutine
subroutine test_dgesvd
    use f77_lapack, only: dgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_dgesvd

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 3
    real(REAL64) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL64) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    real(REAL64), allocatable :: rwork(:)  ! Needed for complex types
    integer :: lwork

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    allocate(work(1))
    lwork = -1  ! Workspace query
        ! For real types, no rwork needed
        call dgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
            call dgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_dgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_dgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_gesvd")

end subroutine
subroutine test_cgesvd
    use f77_lapack, only: cgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_cgesvd

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 3
    complex(REAL32) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL32) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    real(REAL32), allocatable :: rwork(:)  ! Needed for complex types
    integer :: lwork

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    allocate(work(1))
    lwork = -1  ! Workspace query
        ! For complex types, we need rwork as well
        allocate(rwork(5*min(M,N)))
        call cgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, rwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))
            deallocate(rwork)  ! Deallocate workspace from query
            allocate(rwork(5*min(M,N)))  ! Then reallocate for actual call

        A_in = A
            call cgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, rwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_cgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_cgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_gesvd")

end subroutine
subroutine test_zgesvd
    use f77_lapack, only: zgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_zgesvd

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 3
    complex(REAL64) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL64) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    real(REAL64), allocatable :: rwork(:)  ! Needed for complex types
    integer :: lwork

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    allocate(work(1))
    lwork = -1  ! Workspace query
        ! For complex types, we need rwork as well
        allocate(rwork(5*min(M,N)))
        call zgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, rwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))
            deallocate(rwork)  ! Deallocate workspace from query
            allocate(rwork(5*min(M,N)))  ! Then reallocate for actual call

        A_in = A
            call zgesvd('N', 'N', M, N, A_in, M, S, A_in, M, A_in, N, work, lwork, rwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_zgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_zgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf, "different info results for mfi_gesvd")

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
