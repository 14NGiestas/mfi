  program test_hetrf_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: chetrf, zhetrf
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
 call test_chetrf_gpu 
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
        "testing [34mc[0m mfi_hetrf ([32mGPU[0m) against [34mchetrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zhetrf_gpu 
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
        "testing [33mz[0m mfi_hetrf ([32mGPU[0m) against [33mzhetrf[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_chetrf_gpu
    use f77_lapack, only: chetrf, f77_hetrf
    use mfi_blas
    use mfi_lapack, only: mfi_hetrf, mfi_chetrf


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 4
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N)
    complex(REAL32) :: A_original(N,N)  ! Store original matrix for multiple tests
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create a Hermitian matrix for hetrf (complex version)
    ! For simplicity, using a real symmetric matrix (which works for hetrf on real types)
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Store original matrix for multiple tests
    A_original = A

    ! Test f77 interface
    A_in = A_original

    ! Query optimal workspace size
    allocate(work(1))
    lwork = -1
    call chetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N  ! Provide a reasonable default
        lwork = max(1, lwork)
        deallocate(work)
        allocate(work(lwork))

        ! Perform the factorization
        A_in = A_original  ! Reset input matrix to original
        call chetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)

        if (info /= 0) return  ! Skip if factorization failed
    else
        if (allocated(work)) deallocate(work)
        return  ! Skip if workspace query failed
    end if

    ! Test mfi interface (short form) - includes uplo and ipiv parameters
    A_in = A_original  ! Reset to original matrix
    ipiv_in = 0  ! Initialize pivot array
    call mfi_chetrf(A_in, 'U', ipiv_in, info=info_mfi)
    ! Compare info, pivots and factorized matrices for hetrf implementations
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_chetrf")

    ! Test mfi interface (full form) - includes uplo and ipiv parameters
    A_in = A_original  ! Reset to original matrix
    ipiv_in = 0  ! Initialize pivot array
    call mfi_hetrf(A_in, 'U', ipiv_in, info=info_mfi)
    ! Compare info, pivots and factorized matrices for hetrf implementations
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_hetrf")

end subroutine
subroutine test_zhetrf_gpu
    use f77_lapack, only: zhetrf, f77_hetrf
    use mfi_blas
    use mfi_lapack, only: mfi_hetrf, mfi_zhetrf


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 4
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N)
    complex(REAL64) :: A_original(N,N)  ! Store original matrix for multiple tests
    integer :: ipiv_in(N), ipiv_rf(N)
    integer :: info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork

    call mfi_force_gpu()

    ! Create a Hermitian matrix for hetrf (complex version)
    ! For simplicity, using a real symmetric matrix (which works for hetrf on real types)
    A = reshape([2.0_wp, -1.0_wp,  0.0_wp,  0.0_wp, &
                -1.0_wp,  2.0_wp, -1.0_wp,  0.0_wp, &
                 0.0_wp, -1.0_wp,  2.0_wp, -1.0_wp, &
                 0.0_wp,  0.0_wp, -1.0_wp,  2.0_wp], [N,N])

    ! Store original matrix for multiple tests
    A_original = A

    ! Test f77 interface
    A_in = A_original

    ! Query optimal workspace size
    allocate(work(1))
    lwork = -1
    call zhetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = N  ! Provide a reasonable default
        lwork = max(1, lwork)
        deallocate(work)
        allocate(work(lwork))

        ! Perform the factorization
        A_in = A_original  ! Reset input matrix to original
        call zhetrf('U', N, A_in, N, ipiv_in, work, lwork, info)
        A_rf = A_in
        ipiv_rf = ipiv_in
        info_rf = info
        deallocate(work)

        if (info /= 0) return  ! Skip if factorization failed
    else
        if (allocated(work)) deallocate(work)
        return  ! Skip if workspace query failed
    end if

    ! Test mfi interface (short form) - includes uplo and ipiv parameters
    A_in = A_original  ! Reset to original matrix
    ipiv_in = 0  ! Initialize pivot array
    call mfi_zhetrf(A_in, 'U', ipiv_in, info=info_mfi)
    ! Compare info, pivots and factorized matrices for hetrf implementations
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_zhetrf")

    ! Test mfi interface (full form) - includes uplo and ipiv parameters
    A_in = A_original  ! Reset to original matrix
    ipiv_in = 0  ! Initialize pivot array
    call mfi_hetrf(A_in, 'U', ipiv_in, info=info_mfi)
    ! Compare info, pivots and factorized matrices for hetrf implementations
    call assert(info_mfi == info_rf .and. all(abs(A_in - A_rf) < sqrt(epsilon(1.0_wp))) .and. all(ipiv_in == ipiv_rf), &
                "different results for mfi_hetrf")

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


