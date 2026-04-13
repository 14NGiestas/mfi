  program test_trtrs
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: strtrs, dtrtrs, ctrtrs, ztrtrs
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
 call test_strtrs 
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
        "testing [31ms[0m mfi_trtrs ([34mCPU[0m) against [31mstrtrs[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dtrtrs 
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
        "testing [32md[0m mfi_trtrs ([34mCPU[0m) against [32mdtrtrs[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_ctrtrs 
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
        "testing [34mc[0m mfi_trtrs ([34mCPU[0m) against [34mctrtrs[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_ztrtrs 
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
        "testing [33mz[0m mfi_trtrs ([34mCPU[0m) against [33mztrtrs[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_strtrs
    use f77_lapack, only: strtrs, f77_trtrs
    use mfi_blas
    use mfi_lapack, only: mfi_trtrs, mfi_strtrs


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    real(REAL32) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi


    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call strtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_strtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_strtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_strtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_trtrs")

    call mfi_force_cpu()
end subroutine
subroutine test_dtrtrs
    use f77_lapack, only: dtrtrs, f77_trtrs
    use mfi_blas
    use mfi_lapack, only: mfi_trtrs, mfi_dtrtrs


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    real(REAL64) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi


    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call dtrtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_dtrtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_dtrtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dtrtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_trtrs")

    call mfi_force_cpu()
end subroutine
subroutine test_ctrtrs
    use f77_lapack, only: ctrtrs, f77_trtrs
    use mfi_blas
    use mfi_lapack, only: mfi_trtrs, mfi_ctrtrs


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL32) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi


    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call ctrtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_ctrtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_ctrtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ctrtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_trtrs")

    call mfi_force_cpu()
end subroutine
subroutine test_ztrtrs
    use f77_lapack, only: ztrtrs, f77_trtrs
    use mfi_blas
    use mfi_lapack, only: mfi_trtrs, mfi_ztrtrs


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL64) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi


    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call ztrtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_ztrtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_ztrtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ztrtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_trtrs")

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


