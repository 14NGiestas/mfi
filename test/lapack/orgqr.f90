
program test_orgqr
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sorgqr 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orgqr against sorgqr", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dorgqr 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orgqr against dorgqr", t2-t1
end block
contains

subroutine test_sorgqr
    use f77_lapack, only: sorgqr, f77_orgqr
    use mfi_lapack, only: mfi_orgqr, mfi_sorgqr, mfi_geqrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    real(REAL32) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Compute QR factorization first
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return

    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for orgqr
    allocate(work(1))
    lwork = -1
    call sorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call sorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_sorgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sorgqr")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_orgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgqr")

end subroutine
subroutine test_dorgqr
    use f77_lapack, only: dorgqr, f77_orgqr
    use mfi_lapack, only: mfi_orgqr, mfi_dorgqr, mfi_geqrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    real(REAL64) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Compute QR factorization first
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return

    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for orgqr
    allocate(work(1))
    lwork = -1
    call dorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call dorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_dorgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dorgqr")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_orgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgqr")

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
