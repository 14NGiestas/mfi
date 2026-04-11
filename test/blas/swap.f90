 program swap
 use iso_fortran_env
 use mfi_blas
 implicit none
 print '(A)', "testing mfi_swap (CPU) against sswap"
 print '(A)', "testing mfi_swap (CPU) against dswap"
 print '(A)', "testing mfi_swap (CPU) against cswap"
 print '(A)', "testing mfi_swap (CPU) against zswap"
 contains
subroutine test_sswap
    use f77_blas, only: sswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_sswap

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)

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

    x_in = x; y_in = y
    call sswap(N, x_in, 1, y_in, 1)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_sswap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_swap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_dswap
    use f77_blas, only: dswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_dswap

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)

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

    x_in = x; y_in = y
    call dswap(N, x_in, 1, y_in, 1)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_dswap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_swap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_cswap
    use f77_blas, only: cswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_cswap

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)

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

    x_in = x; y_in = y
    call cswap(N, x_in, 1, y_in, 1)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_cswap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_swap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_zswap
    use f77_blas, only: zswap, f77_swap
    use mfi_blas, only: mfi_swap, mfi_zswap

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)

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

    x_in = x; y_in = y
    call zswap(N, x_in, 1, y_in, 1)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_swap(N, x_in, 1, y_in, 1)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_zswap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_swap(x_in, y_in)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:y: mismatch")

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

