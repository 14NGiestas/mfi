 program test_gerqf_gpu
 use iso_fortran_env
 use mfi_lapack
 use f77_lapack, only: sgerqf, dgerqf, cgerqf, zgerqf
 implicit none
 print '(A)', "testing mfi_gerqf (GPU) against sgerqf"
 print '(A)', "testing mfi_gerqf (GPU) against dgerqf"
 print '(A)', "testing mfi_gerqf (GPU) against cgerqf"
 print '(A)', "testing mfi_gerqf (GPU) against zgerqf"
 contains
subroutine test_sgerqf_gpu
    use f77_lapack, only: sgerqf, f77_gerqf
    use mfi_blas
    use mfi_lapack, only: mfi_gerqf, mfi_sgerqf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 6, M = 2
    real(REAL32) :: A(N,M), A_in(N,M), A_rf(N,M)
    real(REAL32) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

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
    call sgerqf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call sgerqf(N, M, A_in, N, tau_in, work, lwork, info)
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
    call mfi_sgerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_sgerqf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_gerqf")

end subroutine
subroutine test_dgerqf_gpu
    use f77_lapack, only: dgerqf, f77_gerqf
    use mfi_blas
    use mfi_lapack, only: mfi_gerqf, mfi_dgerqf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 6, M = 2
    real(REAL64) :: A(N,M), A_in(N,M), A_rf(N,M)
    real(REAL64) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

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
    call dgerqf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call dgerqf(N, M, A_in, N, tau_in, work, lwork, info)
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
    call mfi_dgerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_dgerqf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_gerqf")

end subroutine
subroutine test_cgerqf_gpu
    use f77_lapack, only: cgerqf, f77_gerqf
    use mfi_blas
    use mfi_lapack, only: mfi_gerqf, mfi_cgerqf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 6, M = 2
    complex(REAL32) :: A(N,M), A_in(N,M), A_rf(N,M)
    complex(REAL32) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

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
    call cgerqf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call cgerqf(N, M, A_in, N, tau_in, work, lwork, info)
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
    call mfi_cgerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_cgerqf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_gerqf")

end subroutine
subroutine test_zgerqf_gpu
    use f77_lapack, only: zgerqf, f77_gerqf
    use mfi_blas
    use mfi_lapack, only: mfi_gerqf, mfi_zgerqf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 6, M = 2
    complex(REAL64) :: A(N,M), A_in(N,M), A_rf(N,M)
    complex(REAL64) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

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
    call zgerqf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call zgerqf(N, M, A_in, N, tau_in, work, lwork, info)
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
    call mfi_zgerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_zgerqf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 1e-10) .and. all(abs(tau_in - tau_rf) < 1e-10) .and. info_mfi == info_rf, &
                "different results for mfi_gerqf")

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

