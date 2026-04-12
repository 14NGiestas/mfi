  program test_orgqr
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: sorgqr, dorgqr
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
 call test_sorgqr 
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
        "testing [31ms[0m mfi_orgqr ([34mCPU[0m) against [31msorgqr[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dorgqr 
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
        "testing [32md[0m mfi_orgqr ([34mCPU[0m) against [32mdorgqr[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_sorgqr
    use f77_lapack, only: sorgqr, f77_orgqr
    use mfi_blas
    use mfi_lapack, only: mfi_orgqr, mfi_sorgqr, mfi_geqrf, mfi_gerqf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    real(REAL32) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork


    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Determine which factorization to use based on routine name
    ! Compute QR factorization first for QR routines
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for sorgqr
    allocate(work(1))
    lwork = -1
    call sorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call sorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_sorgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sorgqr")

    ! Test mfi interface (full form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_orgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgqr")

end subroutine
subroutine test_dorgqr
    use f77_lapack, only: dorgqr, f77_orgqr
    use mfi_blas
    use mfi_lapack, only: mfi_orgqr, mfi_dorgqr, mfi_geqrf, mfi_gerqf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    real(REAL64) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork


    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Determine which factorization to use based on routine name
    ! Compute QR factorization first for QR routines
    A_in = A
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for dorgqr
    allocate(work(1))
    lwork = -1
    call dorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call dorgqr(N, N, N, A_in, N, tau_in, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_dorgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dorgqr")

    ! Test mfi interface (full form)
    A_in = A
    ! QR factorization for QR routines
    call mfi_geqrf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_orgqr(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_orgqr")

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


