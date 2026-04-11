

program test_rot
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_srot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against srot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_drot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against drot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_crot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against crot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against zrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_csrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against csrot", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zdrot 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rot against zdrot", t2-t1
end block
contains
subroutine test_srot
    use f77_blas, only: srot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_srot

    integer, parameter :: wp = REAL32
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer, parameter :: N = 20
    real(REAL32) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL32) :: c
    real(REAL32) :: s

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
    call random_number(angle)
    angle = angle * 2.0_wp * pi
    call random_number(x)
    call random_number(y)

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call srot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_srot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_drot
    use f77_blas, only: drot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_drot

    integer, parameter :: wp = REAL64
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer, parameter :: N = 20
    real(REAL64) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL64) :: c
    real(REAL64) :: s

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
    call random_number(angle)
    angle = angle * 2.0_wp * pi
    call random_number(x)
    call random_number(y)

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call drot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_drot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_crot
    use f77_blas, only: crot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_crot

    integer, parameter :: wp = REAL32
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer, parameter :: N = 20
    complex(REAL32) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL32) :: c
    complex(REAL32) :: s

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
    call random_number(angle)
    angle = angle * 2.0_wp * pi
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

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call crot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_crot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_zrot
    use f77_blas, only: zrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_zrot

    integer, parameter :: wp = REAL64
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer, parameter :: N = 20
    complex(REAL64) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL64) :: c
    complex(REAL64) :: s

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
    call random_number(angle)
    angle = angle * 2.0_wp * pi
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

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call zrot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_zrot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_csrot
    use f77_blas, only: csrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_csrot

    integer, parameter :: wp = REAL32
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer, parameter :: N = 20
    complex(REAL32) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL32) :: c
    real(REAL32) :: s

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
    call random_number(angle)
    angle = angle * 2.0_wp * pi
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

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call csrot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_csrot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

end subroutine
subroutine test_zdrot
    use f77_blas, only: zdrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_zdrot

    integer, parameter :: wp = REAL64
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer, parameter :: N = 20
    complex(REAL64) :: x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N)
    real(wp) :: angle
    real(REAL64) :: c
    real(REAL64) :: s

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
    call random_number(angle)
    angle = angle * 2.0_wp * pi
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

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call zdrot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_zdrot(x_in, y_in, c, s)
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
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

