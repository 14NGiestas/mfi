  program test_ormqr
  use iso_fortran_env
  use mfi_lapack
  use f77_lapack, only: sormqr, dormqr
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
 call test_sormqr 
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
        "testing [31ms[0m mfi_ormqr ([34mCPU[0m) against [31msormqr[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dormqr 
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
        "testing [32md[0m mfi_ormqr ([34mCPU[0m) against [32mdormqr[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
subroutine test_sormqr
    use f77_lapack, only: sormqr, f77_ormqr
    use mfi_blas
    use mfi_lapack, only: mfi_ormqr, mfi_sormqr, mfi_geqrf, mfi_gerqf


    integer, parameter :: wp = REAL32
    integer, parameter :: M = 3, N = 2
    real(REAL32) :: A(M,M), A_copy(M,M), C(M,N), C_orig(M,N), C_rf(M,N)
    real(REAL32) :: tau(M)
    integer :: info, info_rf, info_mfi
    real(REAL32), allocatable :: work(:)
    integer :: lwork
    character :: side = 'L', trans = 'N'
    ! For RQ routines, we need temporary matrix for proper factorization
    real(REAL32) :: A_factorize(M,N)  ! Temporary matrix for RQ factorization


    ! Create a test matrix A (will be factorized)
    ! For RQ/QR routines, the factorization matrix A needs appropriate dimensions for multiplication
    ! When using ORMRQ/ORMQR, A should be properly sized for the side operation:
    ! If side='L': A should be KxM where K=min(M,N) and we compute Q*C
    ! If side='R': A should be KxN where K=min(M,N) and we compute C*Q

    ! For the test, we'll start with an appropriately sized square matrix
    ! that can be used for both QR and RQ factorization
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 2.0_wp], [M,M])

    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, &
                 1.0_wp, 1.0_wp], [M,N])

    ! For RQ/QR routines, we need consistent matrix setup for proper testing
    A_copy = A
    ! For RQ routines, we need to be very careful about dimensions
    ! For multiplying C (M×N), the factorization matrix must come from a matrix compatible with the operation
    ! For side='L': A should be from factorizing an M×N matrix
    ! For side='R': A should be from factorizing an M×N matrix
    ! Actually, let me use proper dimensions: if we want to multiply C(M,N) where C is M×N
    ! For side='L': A is result of RQ factorizing a matrix of size (something) × N
    ! For side='R': A is result of RQ factorizing a matrix of size M × (something)

    ! To keep it simple, let's make A the right size for the RQ factorization
    ! Since we're multiplying C(M,N) and A is initially M×M, let me adjust the approach
    ! For RQ factorization of an M×N matrix, we need to work with proper dimensions

    ! Create proper matrices for RQ case - we need to factorize an appropriate matrix
    ! For QR, the original approach is fine
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return

    ! Store original matrices for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface for sormqr - workspace query and then actual call
    allocate(work(1))
    lwork = -1
    call sormqr(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        ! Use the same factorization and original C matrix for the actual computation
        C_rf = C_orig  ! Reset C_rf to original before applying transform
        call sormqr(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form) - use the same factorization and original input
    ! Reset matrices to original values
    C = C_orig
    A_copy = A  ! Get fresh copy to factorize again (which should give same result within numerical tolerance)
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_sormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_sormqr")

    ! Test mfi interface (full form) - use the same factorization and original input
    C = C_orig
    A_copy = A  ! Get fresh copy to factorize again
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_ormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ormqr")

    call mfi_force_cpu()
end subroutine
subroutine test_dormqr
    use f77_lapack, only: dormqr, f77_ormqr
    use mfi_blas
    use mfi_lapack, only: mfi_ormqr, mfi_dormqr, mfi_geqrf, mfi_gerqf


    integer, parameter :: wp = REAL64
    integer, parameter :: M = 3, N = 2
    real(REAL64) :: A(M,M), A_copy(M,M), C(M,N), C_orig(M,N), C_rf(M,N)
    real(REAL64) :: tau(M)
    integer :: info, info_rf, info_mfi
    real(REAL64), allocatable :: work(:)
    integer :: lwork
    character :: side = 'L', trans = 'N'
    ! For RQ routines, we need temporary matrix for proper factorization
    real(REAL64) :: A_factorize(M,N)  ! Temporary matrix for RQ factorization


    ! Create a test matrix A (will be factorized)
    ! For RQ/QR routines, the factorization matrix A needs appropriate dimensions for multiplication
    ! When using ORMRQ/ORMQR, A should be properly sized for the side operation:
    ! If side='L': A should be KxM where K=min(M,N) and we compute Q*C
    ! If side='R': A should be KxN where K=min(M,N) and we compute C*Q

    ! For the test, we'll start with an appropriately sized square matrix
    ! that can be used for both QR and RQ factorization
    A = reshape([3.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 3.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 2.0_wp], [M,M])

    ! Create test matrix C for multiplication
    C = reshape([2.0_wp, 1.0_wp, &
                 1.0_wp, 2.0_wp, &
                 1.0_wp, 1.0_wp], [M,N])

    ! For RQ/QR routines, we need consistent matrix setup for proper testing
    A_copy = A
    ! For RQ routines, we need to be very careful about dimensions
    ! For multiplying C (M×N), the factorization matrix must come from a matrix compatible with the operation
    ! For side='L': A should be from factorizing an M×N matrix
    ! For side='R': A should be from factorizing an M×N matrix
    ! Actually, let me use proper dimensions: if we want to multiply C(M,N) where C is M×N
    ! For side='L': A is result of RQ factorizing a matrix of size (something) × N
    ! For side='R': A is result of RQ factorizing a matrix of size M × (something)

    ! To keep it simple, let's make A the right size for the RQ factorization
    ! Since we're multiplying C(M,N) and A is initially M×M, let me adjust the approach
    ! For RQ factorization of an M×N matrix, we need to work with proper dimensions

    ! Create proper matrices for RQ case - we need to factorize an appropriate matrix
    ! For QR, the original approach is fine
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return

    ! Store original matrices for comparison
    C_orig = C
    C_rf = C

    ! Test f77 interface for dormqr - workspace query and then actual call
    allocate(work(1))
    lwork = -1
    call dormqr(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
    if (info == 0) then
        lwork = int(real(work(1), wp))
        if (lwork <= 0) lwork = M*N
        deallocate(work)
        allocate(work(max(1, lwork)))

        ! Use the same factorization and original C matrix for the actual computation
        C_rf = C_orig  ! Reset C_rf to original before applying transform
        call dormqr(side, trans, M, N, min(M,N), A_copy, M, tau, C_rf, M, work, lwork, info)
        info_rf = info
        deallocate(work)
    else
        if (allocated(work)) deallocate(work)
        return
    end if

    if (info /= 0) return

    ! Test mfi interface (short form) - use the same factorization and original input
    ! Reset matrices to original values
    C = C_orig
    A_copy = A  ! Get fresh copy to factorize again (which should give same result within numerical tolerance)
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_dormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_dormqr")

    ! Test mfi interface (full form) - use the same factorization and original input
    C = C_orig
    A_copy = A  ! Get fresh copy to factorize again
    call mfi_geqrf(A_copy, tau, info=info)
    if (info /= 0) return
    call mfi_ormqr(A_copy, tau, C, side=side, trans=trans, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(C - C_rf) < sqrt(epsilon(1.0_wp))), &
                "different results for mfi_ormqr")

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


