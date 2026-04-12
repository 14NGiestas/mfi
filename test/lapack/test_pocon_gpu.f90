  program test_pocon_gpu
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: spocon, dpocon, cpocon, zpocon
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
 call test_spocon_gpu 
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
        "testing [31ms[0m mfi_pocon ([32mGPU[0m) against [31mspocon[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dpocon_gpu 
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
        "testing [32md[0m mfi_pocon ([32mGPU[0m) against [32mdpocon[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_cpocon_gpu 
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
        "testing [34mc[0m mfi_pocon ([32mGPU[0m) against [34mcpocon[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zpocon_gpu 
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
        "testing [33mz[0m mfi_pocon ([32mGPU[0m) against [33mzpocon[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_spocon_gpu
    use f77_lapack, only: spocon, f77_pocon
    use mfi_blas
    use mfi_lapack, only: mfi_pocon, mfi_spocon


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N)
    real(REAL32) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL32) :: col_sum
    real(REAL32), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    call mfi_force_gpu()

    ! Test f77 interface for pocon
    allocate(work(3*N))  ! Real pocon functions use 3*N for work array
    allocate(iwork(N))
    call spocon('U', N, A, N, anorm, rcond, work, iwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(iwork)) deallocate(iwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_spocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_spocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_pocon")
    end if

end subroutine
subroutine test_dpocon_gpu
    use f77_lapack, only: dpocon, f77_pocon
    use mfi_blas
    use mfi_lapack, only: mfi_pocon, mfi_dpocon


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N)
    real(REAL64) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL64) :: col_sum
    real(REAL64), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    call mfi_force_gpu()

    ! Test f77 interface for pocon
    allocate(work(3*N))  ! Real pocon functions use 3*N for work array
    allocate(iwork(N))
    call dpocon('U', N, A, N, anorm, rcond, work, iwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(iwork)) deallocate(iwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_dpocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_dpocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_pocon")
    end if

end subroutine
subroutine test_cpocon_gpu
    use f77_lapack, only: cpocon, f77_pocon
    use mfi_blas
    use mfi_lapack, only: mfi_pocon, mfi_cpocon


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N)
    real(REAL32) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL32) :: col_sum
    complex(REAL32), allocatable :: work(:)
    real(REAL32), allocatable :: rwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    call mfi_force_gpu()

    ! Test f77 interface for pocon
    allocate(work(2*N))
    allocate(rwork(N))
    call cpocon('U', N, A, N, anorm, rcond, work, rwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(rwork)) deallocate(rwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_cpocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_cpocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_pocon")
    end if

end subroutine
subroutine test_zpocon_gpu
    use f77_lapack, only: zpocon, f77_pocon
    use mfi_blas
    use mfi_lapack, only: mfi_pocon, mfi_zpocon


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N)
    real(REAL64) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL64) :: col_sum
    complex(REAL64), allocatable :: work(:)
    real(REAL64), allocatable :: rwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    call mfi_force_gpu()

    ! Test f77 interface for pocon
    allocate(work(2*N))
    allocate(rwork(N))
    call zpocon('U', N, A, N, anorm, rcond, work, rwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(rwork)) deallocate(rwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_zpocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_zpocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-3, &
                    "different results for mfi_pocon")
    end if

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


