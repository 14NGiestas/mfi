  program rot
  use iso_fortran_env
  use mfi_blas
  implicit none
block
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
 call test_srot 
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        "testing [31ms[0m mfi_rot ([34mCPU[0m) against [31msrot[0m", t_mu, t_ms, t_mn, t_mx, t_n
end block
block
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
 call test_drot 
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        "testing [32md[0m mfi_rot ([34mCPU[0m) against [32mdrot[0m", t_mu, t_ms, t_mn, t_mx, t_n
end block
block
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
 call test_crot 
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        "testing [34mc[0m mfi_rot ([34mCPU[0m) against [34mcrot[0m", t_mu, t_ms, t_mn, t_mx, t_n
end block
block
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
 call test_zrot 
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        "testing [33mz[0m mfi_rot ([34mCPU[0m) against [33mzrot[0m", t_mu, t_ms, t_mn, t_mx, t_n
end block
block
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
 call test_csrot 
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        "testing [33mcs[0m mfi_rot ([34mCPU[0m) against [33mcsrot[0m", t_mu, t_ms, t_mn, t_mx, t_n
end block
block
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
 call test_zdrot 
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        "testing [33mzd[0m mfi_rot ([34mCPU[0m) against [33mzdrot[0m", t_mu, t_ms, t_mn, t_mx, t_n
end block
  contains

  subroutine fmt_time(t, out)
      real, intent(in) :: t
      character(*), intent(out) :: out
      if (t < 1.0e-3) then
          write(out, '(F12.3,"µs")') t * 1.0e6
      else if (t < 1.0) then
          write(out, '(F12.3,"ms")') t * 1.0e3
      else
          write(out, '(F12.3,"s ")') t
      end if
  end subroutine fmt_time

  subroutine test_get_1d(default, n)
      integer, intent(in) :: default
      integer, intent(out) :: n
      integer :: tgs_stat
      character(32) :: tgs_env
      call get_environment_variable('MFI_TEST_ELEMENTS', value=tgs_env, status=tgs_stat)
      if (tgs_stat == 0 .and. len_trim(tgs_env) > 0) then
          read(tgs_env, '(I10)', iostat=tgs_stat) n
      else
          n = default
      end if
      if (n < 1) n = default
  end subroutine test_get_1d

  subroutine test_get_2d(default, n)
      integer, intent(in) :: default
      integer, intent(out) :: n
      integer :: tgs_stat, tgs_val
      character(32) :: tgs_env
      call get_environment_variable('MFI_TEST_ELEMENTS', value=tgs_env, status=tgs_stat)
      if (tgs_stat == 0 .and. len_trim(tgs_env) > 0) then
          read(tgs_env, '(I10)', iostat=tgs_stat) tgs_val
          n = int(sqrt(real(tgs_val)))
      else
          n = int(sqrt(real(default)))
      end if
      if (n < 4) n = 4
  end subroutine test_get_2d
subroutine test_srot
    use f77_blas, only: srot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_srot

    integer, parameter :: wp = REAL32
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: N
    real(REAL32), allocatable :: x(:), y(:), x_in(:), y_in(:), x_rf(:), y_rf(:)
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
    call test_get_1d(2000, N)

    allocate(x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N))
    call random_number(x)
    call random_number(y)

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call srot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), 'f77_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), 'f77_rot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_srot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), 'mfi_srot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), 'mfi_srot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), 'mfi_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), 'mfi_rot:y' // ": mismatch")

    deallocate(x, y, x_in, y_in, x_rf, y_rf)

end subroutine
subroutine test_drot
    use f77_blas, only: drot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_drot

    integer, parameter :: wp = REAL64
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: N
    real(REAL64), allocatable :: x(:), y(:), x_in(:), y_in(:), x_rf(:), y_rf(:)
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
    call test_get_1d(2000, N)

    allocate(x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N))
    call random_number(x)
    call random_number(y)

    c = cos(angle)
    s = sin(angle)

    x_in = x; y_in = y
    call drot(N, x_in, 1, y_in, 1, c, s)
    x_rf = x_in; y_rf = y_in

    x_in = x; y_in = y
    call f77_rot(N, x_in, 1, y_in, 1, c, s)
    call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), 'f77_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), 'f77_rot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_drot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), 'mfi_drot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), 'mfi_drot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), 'mfi_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), 'mfi_rot:y' // ": mismatch")

    deallocate(x, y, x_in, y_in, x_rf, y_rf)

end subroutine
subroutine test_crot
    use f77_blas, only: crot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_crot

    integer, parameter :: wp = REAL32
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: N
    complex(REAL32), allocatable :: x(:), y(:), x_in(:), y_in(:), x_rf(:), y_rf(:)
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
    call test_get_1d(2000, N)

    allocate(x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N))
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
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'f77_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'f77_rot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_crot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_crot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_crot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_rot:y' // ": mismatch")

    deallocate(x, y, x_in, y_in, x_rf, y_rf)

end subroutine
subroutine test_zrot
    use f77_blas, only: zrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_zrot

    integer, parameter :: wp = REAL64
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: N
    complex(REAL64), allocatable :: x(:), y(:), x_in(:), y_in(:), x_rf(:), y_rf(:)
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
    call test_get_1d(2000, N)

    allocate(x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N))
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
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'f77_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'f77_rot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_zrot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_zrot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_zrot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_rot:y' // ": mismatch")

    deallocate(x, y, x_in, y_in, x_rf, y_rf)

end subroutine
subroutine test_csrot
    use f77_blas, only: csrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_csrot

    integer, parameter :: wp = REAL32
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: N
    complex(REAL32), allocatable :: x(:), y(:), x_in(:), y_in(:), x_rf(:), y_rf(:)
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
    call test_get_1d(2000, N)

    allocate(x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N))
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
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'f77_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'f77_rot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_csrot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_csrot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_csrot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_rot:y' // ": mismatch")

    deallocate(x, y, x_in, y_in, x_rf, y_rf)

end subroutine
subroutine test_zdrot
    use f77_blas, only: zdrot, f77_rot
    use mfi_blas, only: mfi_rot, mfi_zdrot

    integer, parameter :: wp = REAL64
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: N
    complex(REAL64), allocatable :: x(:), y(:), x_in(:), y_in(:), x_rf(:), y_rf(:)
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
    call test_get_1d(2000, N)

    allocate(x(N), y(N), x_in(N), y_in(N), x_rf(N), y_rf(N))
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
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'f77_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'f77_rot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_zdrot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_zdrot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_zdrot:y' // ": mismatch")

    x_in = x; y_in = y
    call mfi_rot(x_in, y_in, c, s)
    call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_rot:x' // ": mismatch")
    call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_rot:y' // ": mismatch")

    deallocate(x, y, x_in, y_in, x_rf, y_rf)

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

