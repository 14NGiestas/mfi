
program test_gesvd
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sgesvd 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gesvd against sgesvd", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dgesvd 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gesvd against dgesvd", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cgesvd 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gesvd against cgesvd", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zgesvd 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gesvd against zgesvd", t2-t1
end block
contains

subroutine test_sgesvd
    use f77_lapack, only: sgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_sgesvd

    ! For real types, no rwork needed
    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 3
    real(REAL32) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL32) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
    real(REAL32) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    call sgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call sgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)
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
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_sgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

end subroutine
subroutine test_dgesvd
    use f77_lapack, only: dgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_dgesvd

    ! For real types, no rwork needed
    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 3
    real(REAL64) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL64) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
    real(REAL64) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    call dgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call dgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)
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
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_dgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

end subroutine
subroutine test_cgesvd
    use f77_lapack, only: cgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_cgesvd

    ! For complex types, we need rwork as well
    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 3
    complex(REAL32) :: A(M,N), A_in(M,N), A_rf(M,N)
    complex(REAL32) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
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
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    allocate(rwork(5*min(M,N)))
    call cgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))
            deallocate(rwork)  ! Deallocate workspace from query
            allocate(rwork(5*min(M,N)))  ! Then reallocate for actual call

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call cgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)
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
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_cgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

end subroutine
subroutine test_zgesvd
    use f77_lapack, only: zgesvd, f77_gesvd
    use mfi_lapack, only: mfi_gesvd, mfi_zgesvd

    ! For complex types, we need rwork as well
    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 3
    complex(REAL64) :: A(M,N), A_in(M,N), A_rf(M,N)
    complex(REAL64) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
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
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    allocate(rwork(5*min(M,N)))
    call zgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))
            deallocate(rwork)  ! Deallocate workspace from query
            allocate(rwork(5*min(M,N)))  ! Then reallocate for actual call

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call zgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)
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
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_zgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

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
