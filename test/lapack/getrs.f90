
program test_getrs
    use iso_fortran_env
    implicit none

    call test_sgetrs
    call test_dgetrs
    call test_cgetrs
    call test_zgetrs

contains

subroutine test_sgetrs
    use f77_lapack, only: sgetrs, f77_getrs, f77_getrf, sgetrf
    use mfi_lapack, only: mfi_getrs, mfi_sgetrs

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    real(REAL32) :: A(N,N), A_fact(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N)
    integer :: info, info_rf, info_mfi
    character :: trans = 'N'

    ! Create a test matrix A and factorize it first
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Factorize A first using getrf (so we can solve with getrs)
    A_fact = A
    call sgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getrs
    B_in = B
    call sgetrs(trans, N, NRHS, A_fact, N, ipiv, B_in, N, info)
    B_rf = B_in
    info_rf = info

    ! Test mfi interface (short form)
    A_fact = A  ! Refactorize
    call sgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_sgetrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_sgetrs")

    ! Test mfi interface (full form)
    A_fact = A  ! Refactorize
    call sgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_getrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getrs")

end subroutine
subroutine test_dgetrs
    use f77_lapack, only: dgetrs, f77_getrs, f77_getrf, dgetrf
    use mfi_lapack, only: mfi_getrs, mfi_dgetrs

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    real(REAL64) :: A(N,N), A_fact(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N)
    integer :: info, info_rf, info_mfi
    character :: trans = 'N'

    ! Create a test matrix A and factorize it first
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Factorize A first using getrf (so we can solve with getrs)
    A_fact = A
    call dgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getrs
    B_in = B
    call dgetrs(trans, N, NRHS, A_fact, N, ipiv, B_in, N, info)
    B_rf = B_in
    info_rf = info

    ! Test mfi interface (short form)
    A_fact = A  ! Refactorize
    call dgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_dgetrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_dgetrs")

    ! Test mfi interface (full form)
    A_fact = A  ! Refactorize
    call dgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_getrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getrs")

end subroutine
subroutine test_cgetrs
    use f77_lapack, only: cgetrs, f77_getrs, f77_getrf, cgetrf
    use mfi_lapack, only: mfi_getrs, mfi_cgetrs

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL32) :: A(N,N), A_fact(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N)
    integer :: info, info_rf, info_mfi
    character :: trans = 'N'

    ! Create a test matrix A and factorize it first
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Factorize A first using getrf (so we can solve with getrs)
    A_fact = A
    call cgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getrs
    B_in = B
    call cgetrs(trans, N, NRHS, A_fact, N, ipiv, B_in, N, info)
    B_rf = B_in
    info_rf = info

    ! Test mfi interface (short form)
    A_fact = A  ! Refactorize
    call cgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_cgetrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_cgetrs")

    ! Test mfi interface (full form)
    A_fact = A  ! Refactorize
    call cgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_getrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getrs")

end subroutine
subroutine test_zgetrs
    use f77_lapack, only: zgetrs, f77_getrs, f77_getrf, zgetrf
    use mfi_lapack, only: mfi_getrs, mfi_zgetrs

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL64) :: A(N,N), A_fact(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N)
    integer :: info, info_rf, info_mfi
    character :: trans = 'N'

    ! Create a test matrix A and factorize it first
    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Factorize A first using getrf (so we can solve with getrs)
    A_fact = A
    call zgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getrs
    B_in = B
    call zgetrs(trans, N, NRHS, A_fact, N, ipiv, B_in, N, info)
    B_rf = B_in
    info_rf = info

    ! Test mfi interface (short form)
    A_fact = A  ! Refactorize
    call zgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_zgetrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_zgetrs")

    ! Test mfi interface (full form)
    A_fact = A  ! Refactorize
    call zgetrf(N, N, A_fact, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    B_in = B
    call mfi_getrs(A_fact, ipiv, B_in, info=info_mfi)
    call assert(all(abs(B_in - B_rf) < 1e-8) .and. info_mfi == info_rf, "different results for mfi_getrs")

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
