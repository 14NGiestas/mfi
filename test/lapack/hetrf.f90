
program test_hetrf
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_chetrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_hetrf against chetrf", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zhetrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_hetrf against zhetrf", t2-t1
end block
contains
subroutine test_chetrf
    use f77_lapack, only: chetrf, f77_hetrf
    use mfi_lapack, only: mfi_hetrf, mfi_chetrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 4
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    ! Create a Hermitian matrix for hetrf (complex version)
    ! For simplicity, using a real symmetric matrix (which works for hetrf on real types)
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Test f77 interface
    A_in = A

    ! Query optimal workspace size
    allocate(work(1))
    lwork = -1
    call chetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N  ! Provide a reasonable default
        lwork = max(1, lwork)
        deallocate(work)
        allocate(work(lwork))

        ! Perform the factorization
        A_in = A  ! Reset input matrix
        call chetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)

        if (info /= 0) return  ! Skip if factorization failed
    else
        if (allocated(work)) deallocate(work)
        return  ! Skip if workspace query failed
    end if

    ! Test mfi interface (short form) - includes uplo and ipiv parameters
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_chetrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_chetrf")

    ! Test mfi interface (full form) - includes uplo and ipiv parameters
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_hetrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_hetrf")

end subroutine
subroutine test_zhetrf
    use f77_lapack, only: zhetrf, f77_hetrf
    use mfi_lapack, only: mfi_hetrf, mfi_zhetrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 4
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    ! Create a Hermitian matrix for hetrf (complex version)
    ! For simplicity, using a real symmetric matrix (which works for hetrf on real types)
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Test f77 interface
    A_in = A

    ! Query optimal workspace size
    allocate(work(1))
    lwork = -1
    call zhetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N  ! Provide a reasonable default
        lwork = max(1, lwork)
        deallocate(work)
        allocate(work(lwork))

        ! Perform the factorization
        A_in = A  ! Reset input matrix
        call zhetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)

        if (info /= 0) return  ! Skip if factorization failed
    else
        if (allocated(work)) deallocate(work)
        return  ! Skip if workspace query failed
    end if

    ! Test mfi interface (short form) - includes uplo and ipiv parameters
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_zhetrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_zhetrf")

    ! Test mfi interface (full form) - includes uplo and ipiv parameters
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_hetrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_hetrf")

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
