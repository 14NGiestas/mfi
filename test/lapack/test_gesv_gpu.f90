  program test_gesv_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: sgesv, dgesv, cgesv, zgesv
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
 call test_sgesv_gpu 
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
        "testing [31ms[0m mfi_gesv ([32mGPU[0m) against [31msgesv[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dgesv_gpu 
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
        "testing [32md[0m mfi_gesv ([32mGPU[0m) against [32mdgesv[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_cgesv_gpu 
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
        "testing [34mc[0m mfi_gesv ([32mGPU[0m) against [34mcgesv[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zgesv_gpu 
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
        "testing [33mz[0m mfi_gesv ([32mGPU[0m) against [33mzgesv[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_sgesv_gpu
    use f77_lapack, only: sgesv
    use mfi_blas
    use mfi_lapack, only: mfi_gesv, mfi_sgesv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    real(REAL32) :: A(N,N), A_in(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    call mfi_force_gpu()

    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Test f77 interface for gesv
    A_in = A
    B_in = B
    call sgesv(N, NRHS, A_in, N, ipiv, B_in, N, info)
    B_rf = B_in
    ipiv_rf = ipiv
    info_rf = info

    call assert(info == 0, "f77_sgesv failed")

    ! Test mfi interface (short form)
    A_in = A
    B_in = B
    call mfi_sgesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_sgesv")

    ! Test mfi interface (full form)
    A_in = A
    B_in = B
    call mfi_gesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_gesv")

    call mfi_force_cpu()
end subroutine
subroutine test_dgesv_gpu
    use f77_lapack, only: dgesv
    use mfi_blas
    use mfi_lapack, only: mfi_gesv, mfi_dgesv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    real(REAL64) :: A(N,N), A_in(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    call mfi_force_gpu()

    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Test f77 interface for gesv
    A_in = A
    B_in = B
    call dgesv(N, NRHS, A_in, N, ipiv, B_in, N, info)
    B_rf = B_in
    ipiv_rf = ipiv
    info_rf = info

    call assert(info == 0, "f77_dgesv failed")

    ! Test mfi interface (short form)
    A_in = A
    B_in = B
    call mfi_dgesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_dgesv")

    ! Test mfi interface (full form)
    A_in = A
    B_in = B
    call mfi_gesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_gesv")

    call mfi_force_cpu()
end subroutine
subroutine test_cgesv_gpu
    use f77_lapack, only: cgesv
    use mfi_blas
    use mfi_lapack, only: mfi_gesv, mfi_cgesv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL32) :: A(N,N), A_in(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    call mfi_force_gpu()

    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Test f77 interface for gesv
    A_in = A
    B_in = B
    call cgesv(N, NRHS, A_in, N, ipiv, B_in, N, info)
    B_rf = B_in
    ipiv_rf = ipiv
    info_rf = info

    call assert(info == 0, "f77_cgesv failed")

    ! Test mfi interface (short form)
    A_in = A
    B_in = B
    call mfi_cgesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_cgesv")

    ! Test mfi interface (full form)
    A_in = A
    B_in = B
    call mfi_gesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_gesv")

    call mfi_force_cpu()
end subroutine
subroutine test_zgesv_gpu
    use f77_lapack, only: zgesv
    use mfi_blas
    use mfi_lapack, only: mfi_gesv, mfi_zgesv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL64) :: A(N,N), A_in(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: ipiv(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi

    call mfi_force_gpu()

    A(1,:) = [2.0_wp, 1.0_wp, 1.0_wp]
    A(2,:) = [4.0_wp, 3.0_wp, 3.0_wp]
    A(3,:) = [8.0_wp, 7.0_wp, 9.0_wp]

    B(1,:) = [1.0_wp, 3.0_wp]
    B(2,:) = [1.0_wp, 1.0_wp]
    B(3,:) = [3.0_wp, 1.0_wp]

    ! Test f77 interface for gesv
    A_in = A
    B_in = B
    call zgesv(N, NRHS, A_in, N, ipiv, B_in, N, info)
    B_rf = B_in
    ipiv_rf = ipiv
    info_rf = info

    call assert(info == 0, "f77_zgesv failed")

    ! Test mfi interface (short form)
    A_in = A
    B_in = B
    call mfi_zgesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_zgesv")

    ! Test mfi interface (full form)
    A_in = A
    B_in = B
    call mfi_gesv(A_in, B_in, ipiv, info=info_mfi)
    call assert(info_mfi == info_rf .and. &
                all(abs(B_in - B_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. &
                all(ipiv == ipiv_rf), &
                "different results for mfi_gesv")

    call mfi_force_cpu()
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


