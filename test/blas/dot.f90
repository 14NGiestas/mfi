

program test_dot
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sdot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dot (CPU) against sdot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ddot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_dot (CPU) against ddot", t2-t1
end block
contains
subroutine test_sdot
    use f77_blas, only: sdot, f77_dot
    use mfi_blas, only: mfi_dot, mfi_sdot

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: res, ref
    real(REAL32) :: x(N), y(N)

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
    call random_number(x)
    call random_number(y)

    ref = sdot(N, x, 1, y, 1)

    res = f77_dot(N, x, 1, y, 1)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "f77_dot mismatch")

    res = mfi_sdot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_sdot mismatch")

    res = mfi_dot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_dot mismatch")

end subroutine
subroutine test_ddot
    use f77_blas, only: ddot, f77_dot
    use mfi_blas, only: mfi_dot, mfi_ddot

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: res, ref
    real(REAL64) :: x(N), y(N)

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
    call random_number(x)
    call random_number(y)

    ref = ddot(N, x, 1, y, 1)

    res = f77_dot(N, x, 1, y, 1)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "f77_dot mismatch")

    res = mfi_ddot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_ddot mismatch")

    res = mfi_dot(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_dot mismatch")

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

