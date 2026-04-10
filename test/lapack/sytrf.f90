
program test_sytrf
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_ssytrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_sytrf against ssytrf", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dsytrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_sytrf against dsytrf", t2-t1
end block
contains

subroutine test_ssytrf
    use f77_lapack, only: ssytrf, f77_sytrf
    use mfi_lapack, only: mfi_sytrf, mfi_ssytrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 4
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N), swork(1)
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create a symmetric matrix
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Test f77 interface
    A_in = A
    call ssytrf('U', N, A_in, N, ipiv_in, swork, -1, info)  ! Workspace query
    if (info == 0) then
        lwork = int(swork(1))
        allocate(work(lwork))
        A_in = A
        call ssytrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        return
    end if

    ! Test mfi interface (short form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_ssytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_ssytrf")

    ! Test mfi interface (full form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_sytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_sytrf")

end subroutine
subroutine test_dsytrf
    use f77_lapack, only: dsytrf, f77_sytrf
    use mfi_lapack, only: mfi_sytrf, mfi_dsytrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 4
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N), swork(1)
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create a symmetric matrix
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Test f77 interface
    A_in = A
    call dsytrf('U', N, A_in, N, ipiv_in, swork, -1, info)  ! Workspace query
    if (info == 0) then
        lwork = int(swork(1))
        allocate(work(lwork))
        A_in = A
        call dsytrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        return
    end if

    ! Test mfi interface (short form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_dsytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_dsytrf")

    ! Test mfi interface (full form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_sytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_sytrf")

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
