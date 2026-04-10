
program test_orgrq
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sorgrq 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orgrq against sorgrq", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dorgrq 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_orgrq against dorgrq", t2-t1
end block
contains

subroutine test_sorgrq
    use f77_lapack, only: sorgrq, f77_orgrq
    use mfi_lapack, only: mfi_orgrq, mfi_sorgrq, mfi_geqrf, mfi_gerqf

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

    ! Determine which factorization to use based on routine name
    ! Compute RQ factorization first for RQ routines
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for sorgrq
    allocate(work(1))
    lwork = -1
    call sorgrq(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call sorgrq(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_sorgrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sorgrq")

    ! Test mfi interface (full form)
    A_in = A
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_orgrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgrq")

end subroutine
subroutine test_dorgrq
    use f77_lapack, only: dorgrq, f77_orgrq
    use mfi_lapack, only: mfi_orgrq, mfi_dorgrq, mfi_geqrf, mfi_gerqf

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

    ! Determine which factorization to use based on routine name
    ! Compute RQ factorization first for RQ routines
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for dorgrq
    allocate(work(1))
    lwork = -1
    call dorgrq(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call dorgrq(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_dorgrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dorgrq")

    ! Test mfi interface (full form)
    A_in = A
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_orgrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgrq")

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