  program test_ungrq_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: cungrq, zungrq
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
 call test_cungrq_gpu 
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
        "testing [34mc[0m mfi_ungrq ([32mGPU[0m) against [34mcungrq[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zungrq_gpu 
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
        "testing [33mz[0m mfi_ungrq ([32mGPU[0m) against [33mzungrq[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_cungrq_gpu
    use f77_lapack, only: cungrq, f77_ungrq
    use mfi_blas
    use mfi_lapack, only: mfi_ungrq, mfi_cungrq, mfi_geqrf, mfi_gerqf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    complex(REAL32) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Determine which factorization to use based on routine name
    ! Compute RQ factorization first for RQ routines
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for cungrq
    allocate(work(1))
    lwork = -1
    call cungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call cungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_cungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_cungrq")

    ! Test mfi interface (full form)
    A_in = A
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_ungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungrq")

end subroutine
subroutine test_zungrq_gpu
    use f77_lapack, only: zungrq, f77_ungrq
    use mfi_blas
    use mfi_lapack, only: mfi_ungrq, mfi_zungrq, mfi_geqrf, mfi_gerqf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    complex(REAL64) :: tau_in(N)
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create a test matrix
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 3.0_wp], [N,N])

    ! Determine which factorization to use based on routine name
    ! Compute RQ factorization first for RQ routines
    A_in = A
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    A_rf = A_in  ! Store factorized matrix

    ! Test f77 interface for zungrq
    allocate(work(1))
    lwork = -1
    call zungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A_rf
        call zungrq(N, N, N, A_in, N, tau_in, work, lwork, info)
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
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_zungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_zungrq")

    ! Test mfi interface (full form)
    A_in = A
    ! RQ factorization for RQ routines
    call mfi_gerqf(A_in, tau_in, info=info)
    if (info /= 0) return
    call mfi_ungrq(A_in, tau_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ungrq")

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


