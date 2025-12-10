
program test_getrf
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting getrf tests...'
    call test_sgetrf
    call test_dgetrf
    call test_cgetrf
    call test_zgetrf
    write(*,'(A)') 'All getrf tests completed successfully.'

contains

subroutine test_sgetrf
    use f77_lapack, only: sgetrf, f77_getrf, sgetrf
    use mfi_lapack, only: mfi_getrf, mfi_sgetrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: ipiv(N), ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix for LU factorization
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    ! Test f77 interface for getrf
    A_in = A
    call sgetrf(N, N, A_in, N, ipiv_in, info)
    A_rf = A_in
    ipiv_rf = ipiv_in
    info_rf = info
    
    call assert(info == 0, "f77_sgetrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_sgetrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_sgetrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_getrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_getrf")

end subroutine
subroutine test_dgetrf
    use f77_lapack, only: dgetrf, f77_getrf, dgetrf
    use mfi_lapack, only: mfi_getrf, mfi_dgetrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: ipiv(N), ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix for LU factorization
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    ! Test f77 interface for getrf
    A_in = A
    call dgetrf(N, N, A_in, N, ipiv_in, info)
    A_rf = A_in
    ipiv_rf = ipiv_in
    info_rf = info
    
    call assert(info == 0, "f77_dgetrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_dgetrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_dgetrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_getrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_getrf")

end subroutine
subroutine test_cgetrf
    use f77_lapack, only: cgetrf, f77_getrf, cgetrf
    use mfi_lapack, only: mfi_getrf, mfi_cgetrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: ipiv(N), ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix for LU factorization
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    ! Test f77 interface for getrf
    A_in = A
    call cgetrf(N, N, A_in, N, ipiv_in, info)
    A_rf = A_in
    ipiv_rf = ipiv_in
    info_rf = info
    
    call assert(info == 0, "f77_cgetrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_cgetrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_cgetrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_getrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_getrf")

end subroutine
subroutine test_zgetrf
    use f77_lapack, only: zgetrf, f77_getrf, zgetrf
    use mfi_lapack, only: mfi_getrf, mfi_zgetrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: ipiv(N), ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    ! Create a test matrix for LU factorization
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    ! Test f77 interface for getrf
    A_in = A
    call zgetrf(N, N, A_in, N, ipiv_in, info)
    A_rf = A_in
    ipiv_rf = ipiv_in
    info_rf = info
    
    call assert(info == 0, "f77_zgetrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_zgetrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_zgetrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_getrf(A_in, ipiv_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_getrf")

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
