  program test_potrf
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: spotrf, dpotrf, cpotrf, zpotrf
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
 call test_spotrf 
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
        "testing [31ms[0m mfi_potrf ([34mCPU[0m) against [31mspotrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dpotrf 
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
        "testing [32md[0m mfi_potrf ([34mCPU[0m) against [32mdpotrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_cpotrf 
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
        "testing [34mc[0m mfi_potrf ([34mCPU[0m) against [34mcpotrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zpotrf 
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
        "testing [33mz[0m mfi_potrf ([34mCPU[0m) against [33mzpotrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_spotrf
    use f77_lapack, only: spotrf, f77_potrf, f77_spotrf => spotrf
    use mfi_blas
    use mfi_lapack, only: mfi_potrf, mfi_spotrf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi


    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call spotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_spotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_spotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_spotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_potrf")

    call mfi_force_cpu()
end subroutine
subroutine test_dpotrf
    use f77_lapack, only: dpotrf, f77_potrf, f77_dpotrf => dpotrf
    use mfi_blas
    use mfi_lapack, only: mfi_potrf, mfi_dpotrf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi


    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call dpotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_dpotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_dpotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dpotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_potrf")

    call mfi_force_cpu()
end subroutine
subroutine test_cpotrf
    use f77_lapack, only: cpotrf, f77_potrf, f77_cpotrf => cpotrf
    use mfi_blas
    use mfi_lapack, only: mfi_potrf, mfi_cpotrf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi


    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call cpotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_cpotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_cpotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cpotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_potrf")

    call mfi_force_cpu()
end subroutine
subroutine test_zpotrf
    use f77_lapack, only: zpotrf, f77_potrf, f77_zpotrf => zpotrf
    use mfi_blas
    use mfi_lapack, only: mfi_potrf, mfi_zpotrf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    integer :: info, info_rf, info_mfi


    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Test f77 interface for potrf
    A_in = A
    call zpotrf('U', N, A_in, N, info)
    A_rf = A_in
    info_rf = info

    ! Check if factorization was successful
    call assert(info == 0, "f77_zpotrf failed")

    ! Test mfi interface (short form)
    A_in = A
    call mfi_zpotrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zpotrf")

    ! Test mfi interface (full form)
    A_in = A
    call mfi_potrf(A_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_potrf")

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


