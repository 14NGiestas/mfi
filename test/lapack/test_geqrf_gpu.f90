  program test_geqrf_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: sgeqrf, dgeqrf, cgeqrf, zgeqrf
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
 call test_sgeqrf_gpu 
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
        "testing [31ms[0m mfi_geqrf ([32mGPU[0m) against [31msgeqrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dgeqrf_gpu 
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
        "testing [32md[0m mfi_geqrf ([32mGPU[0m) against [32mdgeqrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_cgeqrf_gpu 
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
        "testing [34mc[0m mfi_geqrf ([32mGPU[0m) against [34mcgeqrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zgeqrf_gpu 
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
        "testing [33mz[0m mfi_geqrf ([32mGPU[0m) against [33mzgeqrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_sgeqrf_gpu
    use f77_lapack, only: sgeqrf, f77_geqrf
    use mfi_blas
    use mfi_lapack, only: mfi_geqrf, mfi_sgeqrf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 6, M = 2
    real(REAL32) :: A(N,M), A_in(N,M), A_rf(N,M)
    real(REAL32) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call sgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call sgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_sgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_sgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_geqrf")

end subroutine
subroutine test_dgeqrf_gpu
    use f77_lapack, only: dgeqrf, f77_geqrf
    use mfi_blas
    use mfi_lapack, only: mfi_geqrf, mfi_dgeqrf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 6, M = 2
    real(REAL64) :: A(N,M), A_in(N,M), A_rf(N,M)
    real(REAL64) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call dgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call dgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_dgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_dgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_geqrf")

end subroutine
subroutine test_cgeqrf_gpu
    use f77_lapack, only: cgeqrf, f77_geqrf
    use mfi_blas
    use mfi_lapack, only: mfi_geqrf, mfi_cgeqrf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 6, M = 2
    complex(REAL32) :: A(N,M), A_in(N,M), A_rf(N,M)
    complex(REAL32) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call cgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call cgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_cgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_cgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_geqrf")

end subroutine
subroutine test_zgeqrf_gpu
    use f77_lapack, only: zgeqrf, f77_geqrf
    use mfi_blas
    use mfi_lapack, only: mfi_geqrf, mfi_zgeqrf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 6, M = 2
    complex(REAL64) :: A(N,M), A_in(N,M), A_rf(N,M)
    complex(REAL64) :: tau_in(min(N,M)), tau_rf(min(N,M))
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    A(1,:) = [  .000000_wp,  2.000000_wp]
    A(2,:) = [ 2.000000_wp, -1.000000_wp]
    A(3,:) = [ 2.000000_wp, -1.000000_wp]
    A(4,:) = [  .000000_wp,  1.500000_wp]
    A(5,:) = [ 2.000000_wp, -1.000000_wp]
    A(6,:) = [ 2.000000_wp, -1.000000_wp]

    ! Test f77 interface for workspace query
    A_in = A
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call zgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        call zgeqrf(N, M, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        tau_rf = tau_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if f77 interface failed

    ! Test mfi interface (short form)
    A_in = A
    call mfi_zgeqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_zgeqrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(tau_in - tau_rf) < sqrt(epsilon(1.0_wp))) .and.&
        & info_mfi == info_rf, &
                "different results for mfi_geqrf")

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


