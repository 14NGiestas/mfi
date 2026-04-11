 program dotu_gpu
 use iso_fortran_env
 use mfi_blas
 implicit none
 print '(A)', "testing mfi_dotu (GPU) against cdotu"
 print '(A)', "testing mfi_dotu (GPU) against zdotu"
 contains
subroutine test_cdotu_gpu
    use f77_blas, only: cdotu, f77_dotu
    use mfi_blas, only: mfi_dotu, mfi_cdotu

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

    ref = cdotu(N, x, 1, y, 1)

    res = f77_dotu(N, x, 1, y, 1)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "f77_dotu mismatch")

    res = mfi_cdotu(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_cdotu mismatch")

    res = mfi_dotu(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_dotu mismatch")

end subroutine
subroutine test_zdotu_gpu
    use f77_blas, only: zdotu, f77_dotu
    use mfi_blas, only: mfi_dotu, mfi_zdotu

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

    ref = zdotu(N, x, 1, y, 1)

    res = f77_dotu(N, x, 1, y, 1)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "f77_dotu mismatch")

    res = mfi_zdotu(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_zdotu mismatch")

    res = mfi_dotu(x, y)
    call assert(abs(ref - res) < sqrt(epsilon(1.0_wp)), "mfi_dotu mismatch")

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

