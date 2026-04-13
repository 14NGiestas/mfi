  program test_sytrf_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: ssytrf, dsytrf
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
 call test_ssytrf_gpu 
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
        "testing [31ms[0m mfi_sytrf ([32mGPU[0m) against [31mssytrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dsytrf_gpu 
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
        "testing [32md[0m mfi_sytrf ([32mGPU[0m) against [32mdsytrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_ssytrf_gpu
    use f77_lapack, only: ssytrf, f77_sytrf
    use mfi_blas
    use mfi_lapack, only: mfi_sytrf, mfi_ssytrf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 4
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N), swork(1)
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create a symmetric matrix
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Test f77 interface
    A_in = A
    call ssytrf('U', N, A_in, N, ipiv_in, swork, -1, info)  ! Workspace query
    if (info == 0) then
        lwork = int(swork(1))
        allocate(work(lwork))
        A_in = A
        call ssytrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        return
    end if

    ! Test mfi interface (short form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_ssytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_ssytrf")

    ! Test mfi interface (full form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_sytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_sytrf")

    call mfi_force_cpu()
end subroutine
subroutine test_dsytrf_gpu
    use f77_lapack, only: dsytrf, f77_sytrf
    use mfi_blas
    use mfi_lapack, only: mfi_sytrf, mfi_dsytrf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 4
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N), swork(1)
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create a symmetric matrix
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Test f77 interface
    A_in = A
    call dsytrf('U', N, A_in, N, ipiv_in, swork, -1, info)  ! Workspace query
    if (info == 0) then
        lwork = int(swork(1))
        allocate(work(lwork))
        A_in = A
        call dsytrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        return
    end if

    ! Test mfi interface (short form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_dsytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_dsytrf")

    ! Test mfi interface (full form)
    A_in = A
    ipiv_in = 0  ! Initialize pivot array
    call mfi_sytrf(A_in, 'U', ipiv_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf) .and. info_mfi == info_rf, "different&
        & results for mfi_sytrf")

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


