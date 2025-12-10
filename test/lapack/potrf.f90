
program test_potrf
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_spotrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_potrf against spotrf", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dpotrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_potrf against dpotrf", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cpotrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_potrf against cpotrf", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zpotrf 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_potrf against zpotrf", t2-t1
end block
contains

subroutine test_spotrf
    use f77_lapack, only: spotrf, f77_potrf, f77_spotrf => spotrf
    use mfi_lapack, only: mfi_potrf, mfi_spotrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call spotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_spotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_spotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_spotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_potrf")

end subroutine
subroutine test_dpotrf
    use f77_lapack, only: dpotrf, f77_potrf, f77_dpotrf => dpotrf
    use mfi_lapack, only: mfi_potrf, mfi_dpotrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call dpotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_dpotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_dpotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_dpotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_potrf")

end subroutine
subroutine test_cpotrf
    use f77_lapack, only: cpotrf, f77_potrf, f77_cpotrf => cpotrf
    use mfi_lapack, only: mfi_potrf, mfi_cpotrf

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call cpotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_cpotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_cpotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_cpotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_potrf")

end subroutine
subroutine test_zpotrf
    use f77_lapack, only: zpotrf, f77_potrf, f77_zpotrf => zpotrf
    use mfi_lapack, only: mfi_potrf, mfi_zpotrf

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call zpotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_zpotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_zpotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_zpotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < 1e-10), &
                "different results for mfi_potrf")

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
