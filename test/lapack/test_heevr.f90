  program test_heevr
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: cheevr, zheevr
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
 call test_cheevr 
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
        "testing [34mc[0m mfi_heevr ([34mCPU[0m) against [34mcheevr[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zheevr 
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
        "testing [33mz[0m mfi_heevr ([34mCPU[0m) against [33mzheevr[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_cheevr
    use f77_lapack, only: cheevr, f77_heevr
    use mfi_blas
    use mfi_lapack, only: mfi_heevr, mfi_cheevr


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_copy(N,N), Z(N,N), Z_rf(N,N)
    real(REAL32) :: W(N), W_rf(N)
    complex(REAL32), allocatable :: work(:)
    real(REAL32), allocatable :: rwork(:)
    integer, allocatable :: iwork(:), isuppz(:)
    integer :: info, info_rf, info_mfi, M, M_rf
    character :: jobz = 'V', uplo = 'U', range = 'A'
    real(REAL32) :: vl = 0.0_wp, vu = 0.0_wp
    integer :: il = 1, iu = N
    real(REAL32) :: abstol = 0.0_wp

    ! Create a Hermitian test matrix
    A = 0.0_wp
    A(1,1) = cmplx(3.0_wp, 0.0_wp); A(1,2) = cmplx(1.0_wp, 0.5_wp); A(1,3) = cmplx(1.0_wp, -0.3_wp)
    A(2,1) = cmplx(1.0_wp, -0.5_wp); A(2,2) = cmplx(2.0_wp, 0.0_wp); A(2,3) = cmplx(0.5_wp, 0.2_wp)
    A(3,1) = cmplx(1.0_wp, 0.3_wp); A(3,2) = cmplx(0.5_wp, -0.2_wp); A(3,3) = cmplx(4.0_wp, 0.0_wp)
    
    A_copy = A


    ! Test f77 interface
    allocate(work(26*N), rwork(24*N), iwork(10*N), isuppz(N*2))
    call cheevr(jobz, range, uplo, N, A_copy, N, vl, vu, il, iu, abstol, &
                 M_rf, W_rf, Z_rf, N, isuppz, work, size(work), rwork, size(rwork), iwork, size(iwork), info_rf)

    if (info_rf /= 0) then
        deallocate(work, rwork, iwork, isuppz)
        return
    end if

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    call mfi_cheevr(A_copy, W, jobz=jobz, uplo=uplo, range=range, &
                     m=M, z=Z, isuppz=isuppz, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(W - W_rf) < sqrt(epsilon(1.0_wp))), &
                "different eigenvalues for mfi_cheevr")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    call mfi_heevr(A_copy, W, jobz=jobz, uplo=uplo, range=range, &
                 m=M, z=Z, isuppz=isuppz, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(W - W_rf) < sqrt(epsilon(1.0_wp))), &
                "different eigenvalues for mfi_heevr")

    ! Clean up
    deallocate(work, rwork, iwork, isuppz)

    call mfi_force_cpu()
end subroutine
subroutine test_zheevr
    use f77_lapack, only: zheevr, f77_heevr
    use mfi_blas
    use mfi_lapack, only: mfi_heevr, mfi_zheevr


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_copy(N,N), Z(N,N), Z_rf(N,N)
    real(REAL64) :: W(N), W_rf(N)
    complex(REAL64), allocatable :: work(:)
    real(REAL64), allocatable :: rwork(:)
    integer, allocatable :: iwork(:), isuppz(:)
    integer :: info, info_rf, info_mfi, M, M_rf
    character :: jobz = 'V', uplo = 'U', range = 'A'
    real(REAL64) :: vl = 0.0_wp, vu = 0.0_wp
    integer :: il = 1, iu = N
    real(REAL64) :: abstol = 0.0_wp

    ! Create a Hermitian test matrix
    A = 0.0_wp
    A(1,1) = cmplx(3.0_wp, 0.0_wp); A(1,2) = cmplx(1.0_wp, 0.5_wp); A(1,3) = cmplx(1.0_wp, -0.3_wp)
    A(2,1) = cmplx(1.0_wp, -0.5_wp); A(2,2) = cmplx(2.0_wp, 0.0_wp); A(2,3) = cmplx(0.5_wp, 0.2_wp)
    A(3,1) = cmplx(1.0_wp, 0.3_wp); A(3,2) = cmplx(0.5_wp, -0.2_wp); A(3,3) = cmplx(4.0_wp, 0.0_wp)
    
    A_copy = A


    ! Test f77 interface
    allocate(work(26*N), rwork(24*N), iwork(10*N), isuppz(N*2))
    call zheevr(jobz, range, uplo, N, A_copy, N, vl, vu, il, iu, abstol, &
                 M_rf, W_rf, Z_rf, N, isuppz, work, size(work), rwork, size(rwork), iwork, size(iwork), info_rf)

    if (info_rf /= 0) then
        deallocate(work, rwork, iwork, isuppz)
        return
    end if

    ! Test mfi interface (short form)
    A_copy = A  ! Reset to original
    call mfi_zheevr(A_copy, W, jobz=jobz, uplo=uplo, range=range, &
                     m=M, z=Z, isuppz=isuppz, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(W - W_rf) < sqrt(epsilon(1.0_wp))), &
                "different eigenvalues for mfi_zheevr")

    ! Test mfi interface (full form)
    A_copy = A  ! Reset to original
    call mfi_heevr(A_copy, W, jobz=jobz, uplo=uplo, range=range, &
                 m=M, z=Z, isuppz=isuppz, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(W - W_rf) < sqrt(epsilon(1.0_wp))), &
                "different eigenvalues for mfi_heevr")

    ! Clean up
    deallocate(work, rwork, iwork, isuppz)

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


