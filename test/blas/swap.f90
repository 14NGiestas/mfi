



program test_swap
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against sswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against dswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against cswap", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zswap 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_swap against zswap", t2-t1
end block
contains
subroutine test_sswap
    use f77_blas, only: sswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_sswap

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(REAL32) :: rnd(N)

    real(REAL32) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(X)
    call random_number(Y)

    x_in = x
    y_in = y
    ! The test is always against the original
    call sswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_sswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_dswap
    use f77_blas, only: dswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_dswap

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(REAL64) :: rnd(N)

    real(REAL64) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(X)
    call random_number(Y)

    x_in = x
    y_in = y
    ! The test is always against the original
    call dswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_dswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_cswap
    use f77_blas, only: cswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_cswap

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20

    real(REAL32) :: rnd(N)

    complex(REAL32) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    x_in = x
    y_in = y
    ! The test is always against the original
    call cswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_cswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

end subroutine
subroutine test_zswap
    use f77_blas, only: zswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_zswap

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20

    real(REAL64) :: rnd(N)

    complex(REAL64) :: x(N),    y(N),    &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)

    call random_number(rnd)
    x%re = rnd
    call random_number(rnd)
    x%im = rnd
    call random_number(rnd)
    y%re = rnd
    call random_number(rnd)
    y%im = rnd

    x_in = x
    y_in = y
    ! The test is always against the original
    call zswap(N, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = x
    y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_zswap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

    x_in = x
    y_in = y
    call mfi_swap(x_in, y_in)
    call assert(all(x_in == x_rf) .and. all(y_in == y_rf), "different results")

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

