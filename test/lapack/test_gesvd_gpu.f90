  program test_gesvd_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: sgesvd, dgesvd, cgesvd, zgesvd
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
 call test_sgesvd_gpu 
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
        "testing [31ms[0m mfi_gesvd ([32mGPU[0m) against [31msgesvd[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dgesvd_gpu 
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
        "testing [32md[0m mfi_gesvd ([32mGPU[0m) against [32mdgesvd[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_cgesvd_gpu 
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
        "testing [34mc[0m mfi_gesvd ([32mGPU[0m) against [34mcgesvd[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zgesvd_gpu 
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
        "testing [33mz[0m mfi_gesvd ([32mGPU[0m) against [33mzgesvd[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_sgesvd_gpu
    use f77_lapack, only: sgesvd, f77_gesvd
    use mfi_blas
    use mfi_lapack, only: mfi_gesvd, mfi_sgesvd

    ! For real types, no rwork needed
    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 3
    real(REAL32) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL32) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
    real(REAL32) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    call sgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call sgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_sgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_sgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

    call mfi_force_cpu()
end subroutine
subroutine test_dgesvd_gpu
    use f77_lapack, only: dgesvd, f77_gesvd
    use mfi_blas
    use mfi_lapack, only: mfi_gesvd, mfi_dgesvd

    ! For real types, no rwork needed
    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 3
    real(REAL64) :: A(M,N), A_in(M,N), A_rf(M,N)
    real(REAL64) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
    real(REAL64) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    call dgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call dgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_dgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_dgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

    call mfi_force_cpu()
end subroutine
subroutine test_cgesvd_gpu
    use f77_lapack, only: cgesvd, f77_gesvd
    use mfi_blas
    use mfi_lapack, only: mfi_gesvd, mfi_cgesvd

    ! For complex types, we need rwork as well
    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 3
    complex(REAL32) :: A(M,N), A_in(M,N), A_rf(M,N)
    complex(REAL32) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
    real(REAL32) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    real(REAL32), allocatable :: rwork(:)  ! Needed for complex types
    integer :: lwork

    call mfi_force_gpu()

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    allocate(rwork(5*min(M,N)))
    call cgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))
            deallocate(rwork)  ! Deallocate workspace from query
            allocate(rwork(5*min(M,N)))  ! Then reallocate for actual call

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call cgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_cgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_cgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

    call mfi_force_cpu()
end subroutine
subroutine test_zgesvd_gpu
    use f77_lapack, only: zgesvd, f77_gesvd
    use mfi_blas
    use mfi_lapack, only: mfi_gesvd, mfi_zgesvd

    ! For complex types, we need rwork as well
    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 3
    complex(REAL64) :: A(M,N), A_in(M,N), A_rf(M,N)
    complex(REAL64) :: U_temp(M,M), VT_temp(N,N)  ! Temporary arrays to avoid intent conflicts
    real(REAL64) :: S(min(M,N)), S_rf(min(M,N)), S_mfi(min(M,N))
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    real(REAL64), allocatable :: rwork(:)  ! Needed for complex types
    integer :: lwork

    call mfi_force_gpu()

    ! Create test matrix
    A(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    A(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    A(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]

    ! Test f77 interface (just get S values, not U/V)
    A_in = A
    U_temp = 0.0_wp
    VT_temp = 0.0_wp
    allocate(work(1))
    lwork = -1  ! Workspace query
    allocate(rwork(5*min(M,N)))
    call zgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)

    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))
            deallocate(rwork)  ! Deallocate workspace from query
            allocate(rwork(5*min(M,N)))  ! Then reallocate for actual call

        A_in = A
        U_temp = 0.0_wp  ! Initialize U_temp
        VT_temp = 0.0_wp  ! Initialize VT_temp
            call zgesvd('N', 'N', M, N, A_in, M, S, U_temp, M, VT_temp, N, work, lwork, rwork, info)
        A_rf = A_in
        S_rf = S
        info_rf = info
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
    else
        ! If workspace query failed, skip this test
        deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
        return
    end if

    ! Test mfi interface
    A_in = A
    call mfi_zgesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_zgesvd")

    A_in = A
    call mfi_gesvd(A_in, S_mfi, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(S_mfi - S_rf) < sqrt(epsilon(1.0_wp))), "different results for mfi_gesvd")

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


