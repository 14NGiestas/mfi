

program test_nrm2
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_snrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 (CPU) against snrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dnrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 (CPU) against dnrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_scnrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 (CPU) against scnrm2", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dznrm2 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_nrm2 (CPU) against dznrm2", t2-t1
end block
contains
subroutine test_snrm2
    use f77_blas, only: snrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_snrm2

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = snrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_snrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

    call random_number(array)
    res(1) = snrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_snrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

end subroutine
subroutine test_dnrm2
    use f77_blas, only: dnrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_dnrm2

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = dnrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_dnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

    call random_number(array)
    res(1) = dnrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_dnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

end subroutine
subroutine test_scnrm2
    use f77_blas, only: scnrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_scnrm2

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = scnrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_scnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im, kind=REAL32)
end block
    res(1) = scnrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_scnrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

end subroutine
subroutine test_dznrm2
    use f77_blas, only: dznrm2, f77_nrm2
    use mfi_blas, only: mfi_nrm2, mfi_dznrm2

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = dznrm2(N, array, 1)
    res(3) = f77_nrm2(N, array, 1)
    res(2) = mfi_dznrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im, kind=REAL64)
end block
    res(1) = dznrm2(N, array, 1)
    res(2) = f77_nrm2(N, array, 1)
    res(3) = mfi_dznrm2(array)
    res(4) = mfi_nrm2(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

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

