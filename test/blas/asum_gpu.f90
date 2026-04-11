 program asum_gpu
 use iso_fortran_env
 use mfi_blas
 implicit none
 print '(A)', "testing mfi_asum (GPU) against sasum"
 print '(A)', "testing mfi_asum (GPU) against dasum"
 print '(A)', "testing mfi_asum (GPU) against scasum"
 print '(A)', "testing mfi_asum (GPU) against dzasum"
 contains
subroutine test_sasum_gpu
    use f77_blas, only: sasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_sasum

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: ii

    call mfi_force_gpu()

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
    res(1) = sasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_sasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

    call random_number(array)
    res(1) = sasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_sasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

end subroutine
subroutine test_dasum_gpu
    use f77_blas, only: dasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_dasum

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: ii

    call mfi_force_gpu()

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
    res(1) = dasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_dasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

    call random_number(array)
    res(1) = dasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_dasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

end subroutine
subroutine test_scasum_gpu
    use f77_blas, only: scasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_scasum

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: array(N)
    real(REAL32) :: res(4)
    integer :: ii

    call mfi_force_gpu()

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
    res(1) = scasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_scasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im, kind=REAL32)
end block
    res(1) = scasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_scasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "random array mismatch")

end subroutine
subroutine test_dzasum_gpu
    use f77_blas, only: dzasum, f77_asum
    use mfi_blas, only: mfi_asum, mfi_dzasum

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: array(N)
    real(REAL64) :: res(4)
    integer :: ii

    call mfi_force_gpu()

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
    res(1) = dzasum(N, array, 1)
    res(3) = f77_asum(N, array, 1)
    res(2) = mfi_dzasum(array)
    res(4) = mfi_asum(array)
    call assert(all(abs(res - res(1)) < sqrt(epsilon(1.0_wp))), "sequential array mismatch")

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im, kind=REAL64)
end block
    res(1) = dzasum(N, array, 1)
    res(2) = f77_asum(N, array, 1)
    res(3) = mfi_dzasum(array)
    res(4) = mfi_asum(array)
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

