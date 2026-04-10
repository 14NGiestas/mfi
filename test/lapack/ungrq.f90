
program test_ungrq
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cungrq 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ungrq against cungrq", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zungrq 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_ungrq against zungrq", t2-t1
end block
contains

subroutine test_cungrq
    use f77_lapack, only: cungrq, f77_ungrq
    use mfi_lapack, only: mfi_ungrq, mfi_cungrq, mfi_geqrf, mfi_gerqf

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
    ! Compute RQ factorization first for RQ routines
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for cungrq
    allocate(work(1))
    lwork = -1
    call cungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call cungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    call mfi_cungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cungrq")

    ! Test mfi interface (full form)
    A_in = A
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_ungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungrq")

end subroutine
subroutine test_zungrq
    use f77_lapack, only: zungrq, f77_ungrq
    use mfi_lapack, only: mfi_ungrq, mfi_zungrq, mfi_geqrf, mfi_gerqf

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
    ! Compute RQ factorization first for RQ routines
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for zungrq
    allocate(work(1))
    lwork = -1
    call zungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call zungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    call mfi_zungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zungrq")

    ! Test mfi interface (full form)
    A_in = A
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_ungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungrq")

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