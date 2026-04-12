  program test_getri
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: sgetri, dgetri, cgetri, zgetri
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
 call test_sgetri 
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
        "testing [31ms[0m mfi_getri ([34mCPU[0m) against [31msgetri[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dgetri 
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
        "testing [32md[0m mfi_getri ([34mCPU[0m) against [32mdgetri[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_cgetri 
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
        "testing [34mc[0m mfi_getri ([34mCPU[0m) against [34mcgetri[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_zgetri 
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
        "testing [33mz[0m mfi_getri ([34mCPU[0m) against [33mzgetri[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_sgetri
    use f77_lapack, only: sgetri, f77_getri, f77_getrf, sgetrf
    use mfi_blas
    use mfi_lapack, only: mfi_getri, mfi_sgetri


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork


    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call sgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call sgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call sgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_sgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_sgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_getri")

end subroutine
subroutine test_dgetri
    use f77_lapack, only: dgetri, f77_getri, f77_getrf, dgetrf
    use mfi_blas
    use mfi_lapack, only: mfi_getri, mfi_dgetri


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork


    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call dgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call dgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call dgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_dgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_dgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_getri")

end subroutine
subroutine test_cgetri
    use f77_lapack, only: cgetri, f77_getri, f77_getrf, cgetrf
    use mfi_blas
    use mfi_lapack, only: mfi_getri, mfi_cgetri


    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    complex(REAL32), allocatable :: work(:)
    integer :: lwork


    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call cgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call cgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call cgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_cgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_cgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_getri")

end subroutine
subroutine test_zgetri
    use f77_lapack, only: zgetri, f77_getri, f77_getrf, zgetrf
    use mfi_blas
    use mfi_lapack, only: mfi_getri, mfi_zgetri


    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N), A_in(N,N), A_rf(N,N), A_copy(N,N)
    integer :: ipiv(N), info, info_rf, info_mfi
    complex(REAL64), allocatable :: work(:)
    integer :: lwork


    ! Create an invertible matrix
    A_copy(1,:) = [2.0_wp, 0.0_wp, 1.0_wp]
    A_copy(2,:) = [1.0_wp, 2.0_wp, 0.0_wp]
    A_copy(3,:) = [0.0_wp, 1.0_wp, 2.0_wp]

    ! Factor the matrix first with getrf to prepare for getri
    A = A_copy
    call zgetrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip test if factorization failed

    ! Test f77 interface for getri
    A_in = A  ! A is now the factorized matrix
    allocate(work(1))  ! Small workspace for query
    lwork = -1  ! Workspace query
    call zgetri(N, A_in, N, ipiv, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))  ! Get workspace size
        deallocate(work)
        allocate(work(max(1, lwork)))

        A_in = A  ! Use the factored matrix again
        call zgetri(N, A_in, N, ipiv, work, lwork, info)
        A_rf = A_in
        info_rf = info
        deallocate(work)
    else
        return ! Skip if workspace query failed
    end if

    if (info /= 0) return  ! Skip test if getri failed

    ! Test mfi interface (short form)
    ! We need to factorize again as getri overwrites the input
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_zgetri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_zgetri")

    ! Test mfi interface (full form)
    A = A_copy
    call f77_getrf(N, N, A, N, ipiv, info)
    if (info /= 0) return  ! Skip if factorization failed again

    A_in = A
    call mfi_getri(A_in, ipiv, info=info_mfi)
    call assert(all(abs(A_in - A_rf) < 10.0 * sqrt(epsilon(1.0_wp))) .and. info_mfi == info_rf, "different results for mfi_getri")

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


