
program test_ungqr
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cungqr 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ungqr against cungqr", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zungqr 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ungqr against zungqr", t2-t1
end block
contains

subroutine test_cungqr
    use f77_lapack, only: cungqr, f77_ungqr
    use mfi_lapack, only: mfi_ungqr, mfi_cungqr, mfi_geqrf, mfi_gerqf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    complex(REAL32) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Determine which factorization to use based on routine name
    ! Compute QR factorization first for QR routines
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for cungqr
    allocate(work(1))
    lwork = -1
    call cungqr(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call cungqr(N, N, N, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_cungqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cungqr")

    ! Test mfi interface (full form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_ungqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungqr")

end subroutine
subroutine test_zungqr
    use f77_lapack, only: zungqr, f77_ungqr
    use mfi_lapack, only: mfi_ungqr, mfi_zungqr, mfi_geqrf, mfi_gerqf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    complex(REAL64) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Determine which factorization to use based on routine name
    ! Compute QR factorization first for QR routines
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for zungqr
    allocate(work(1))
    lwork = -1
    call zungqr(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call zungqr(N, N, N, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_zungqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zungqr")

    ! Test mfi interface (full form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_ungqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungqr")

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