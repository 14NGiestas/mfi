

program test_dotc
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_cdotc 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dotc against cdotc", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdotc 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dotc against zdotc", t2-t1
end block
contains
subroutine test_cdotc
    use f77_blas, only: cdotc, f77_dotc
    use mfi_blas, only: mfi_dotc, mfi_cdotc

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: res, ref
    complex(REAL32) :: x(N), y(N)

block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    x = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im, kind=REAL32)
end block

    ref = cdotc(N, x, 1, y, 1)

    res = f77_dotc(N, x, 1, y, 1)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "f77_dotc mismatch")

    res = mfi_cdotc(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_cdotc mismatch")

    res = mfi_dotc(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_dotc mismatch")

end subroutine
subroutine test_zdotc
    use f77_blas, only: zdotc, f77_dotc
    use mfi_blas, only: mfi_dotc, mfi_zdotc

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: res, ref
    complex(REAL64) :: x(N), y(N)

block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    x = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im, kind=REAL64)
end block

    ref = zdotc(N, x, 1, y, 1)

    res = f77_dotc(N, x, 1, y, 1)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "f77_dotc mismatch")

    res = mfi_zdotc(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_zdotc mismatch")

    res = mfi_dotc(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_dotc mismatch")

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

