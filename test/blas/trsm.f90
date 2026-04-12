  program trsm
  use iso_fortran_env
  use mfi_blas
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
 call test_strsm 
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
        "testing [31ms[0m mfi_trsm ([34mCPU[0m) against [31mstrsm[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_dtrsm 
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
        "testing [32md[0m mfi_trsm ([34mCPU[0m) against [32mdtrsm[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_ctrsm 
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
        "testing [34mc[0m mfi_trsm ([34mCPU[0m) against [34mctrsm[0m", t_mu, t_ms, t_mn, t_mx, t_n
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
 call test_ztrsm 
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
        "testing [33mz[0m mfi_trsm ([34mCPU[0m) against [33mztrsm[0m", t_mu, t_ms, t_mn, t_mx, t_n
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

  subroutine test_get_1d(default, n)
      integer, intent(in) :: default
      integer, intent(out) :: n
      integer :: tgs_stat
      character(32) :: tgs_env
      call get_environment_variable('MFI_TEST_ELEMENTS', value=tgs_env, status=tgs_stat)
      if (tgs_stat == 0 .and. len_trim(tgs_env) > 0) then
          read(tgs_env, '(I10)', iostat=tgs_stat) n
      else
          n = default
      end if
      if (n < 1) n = default
  end subroutine test_get_1d

  subroutine test_get_2d(default, n)
      integer, intent(in) :: default
      integer, intent(out) :: n
      integer :: tgs_stat, tgs_val
      character(32) :: tgs_env
      call get_environment_variable('MFI_TEST_ELEMENTS', value=tgs_env, status=tgs_stat)
      if (tgs_stat == 0 .and. len_trim(tgs_env) > 0) then
          read(tgs_env, '(I10)', iostat=tgs_stat) tgs_val
          n = int(sqrt(real(tgs_val)))
      else
          n = int(sqrt(real(default)))
      end if
      if (n < 4) n = 4
  end subroutine test_get_2d
subroutine test_strsm
    use f77_blas, only: strsm, f77_trsm
    use mfi_blas, only: mfi_trsm, mfi_force_gpu, mfi_force_cpu

    integer, parameter :: wp = REAL32
    integer :: M, N
    real(REAL32), allocatable :: A(:,:), B(:,:), B_in(:,:), B_rf(:,:)
    real(REAL32) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag
    integer :: ir, ic

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
    call test_get_2d(2000, M)
    N = M

    allocate(A(M,M), B(M,N), B_in(M,N), B_rf(M,N))
    ! Generate a well-conditioned triangular matrix
    call random_number(A)
    call random_number(B)
    call random_number(alpha)
    ! Zero out half, scale off-diagonal, set diagonal to [1,2]
    do ir = 1, M
        do ic = 1, M
            if (ic > ir) then
                A(ir, ic) = 0.0_wp
            else if (ic == ir) then
                A(ir, ic) = 1.0_wp + real(ir, wp) * 0.1_wp
            else
                A(ir, ic) = (A(ir, ic) - 0.5_wp) * 0.2_wp
            end if
        end do
    end do

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side   = sides(i_side)
        uplo   = uplos(i_uplo)
        transa = transas(i_transa)
        diag   = diags(i_diag)

        B_rf = B
        call strsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
    call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL32)), 'mfi_trsm:' // side // uplo // transa // diag // ": mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
    call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL32)), 'GPU:mfi_trsm:' // side // uplo // transa // diag // ":&
        & mismatch")
#endif

    end do
    end do
    end do
    end do
    deallocate(A, B, B_in, B_rf)
end subroutine
subroutine test_dtrsm
    use f77_blas, only: dtrsm, f77_trsm
    use mfi_blas, only: mfi_trsm, mfi_force_gpu, mfi_force_cpu

    integer, parameter :: wp = REAL64
    integer :: M, N
    real(REAL64), allocatable :: A(:,:), B(:,:), B_in(:,:), B_rf(:,:)
    real(REAL64) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag
    integer :: ir, ic

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
    call test_get_2d(2000, M)
    N = M

    allocate(A(M,M), B(M,N), B_in(M,N), B_rf(M,N))
    ! Generate a well-conditioned triangular matrix
    call random_number(A)
    call random_number(B)
    call random_number(alpha)
    ! Zero out half, scale off-diagonal, set diagonal to [1,2]
    do ir = 1, M
        do ic = 1, M
            if (ic > ir) then
                A(ir, ic) = 0.0_wp
            else if (ic == ir) then
                A(ir, ic) = 1.0_wp + real(ir, wp) * 0.1_wp
            else
                A(ir, ic) = (A(ir, ic) - 0.5_wp) * 0.2_wp
            end if
        end do
    end do

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side   = sides(i_side)
        uplo   = uplos(i_uplo)
        transa = transas(i_transa)
        diag   = diags(i_diag)

        B_rf = B
        call dtrsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
    call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL64)), 'mfi_trsm:' // side // uplo // transa // diag // ": mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
    call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL64)), 'GPU:mfi_trsm:' // side // uplo // transa // diag // ":&
        & mismatch")
#endif

    end do
    end do
    end do
    end do
    deallocate(A, B, B_in, B_rf)
end subroutine
subroutine test_ctrsm
    use f77_blas, only: ctrsm, f77_trsm
    use mfi_blas, only: mfi_trsm, mfi_force_gpu, mfi_force_cpu

    integer, parameter :: wp = REAL32
    integer :: M, N
    complex(REAL32), allocatable :: A(:,:), B(:,:), B_in(:,:), B_rf(:,:)
    complex(REAL32) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag
    integer :: ir, ic

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
    call test_get_2d(2000, M)
    N = M

    allocate(A(M,M), B(M,N), B_in(M,N), B_rf(M,N))
    ! Generate a well-conditioned triangular matrix
block
    real(REAL32) :: re(M,M)
    real(REAL32) :: im(M,M)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re(M,N)
    real(REAL32) :: im(M,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL32)
end block
    ! Zero out half, scale off-diagonal, set diagonal to [1,2]
    do ir = 1, M
        do ic = 1, M
            if (ic > ir) then
                A(ir, ic) = 0.0_wp
            else if (ic == ir) then
                A(ir, ic) = 1.0_wp + real(ir, wp) * 0.1_wp
            else
                A(ir, ic) = (A(ir, ic) - 0.5_wp) * 0.2_wp
            end if
        end do
    end do

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side   = sides(i_side)
        uplo   = uplos(i_uplo)
        transa = transas(i_transa)
        diag   = diags(i_diag)

        B_rf = B
        call ctrsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
    call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'mfi_trsm:' // side // uplo // transa // diag // ":&
        & mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
    call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), 'GPU:mfi_trsm:' // side // uplo // transa // diag //&
        & ": mismatch")
#endif

    end do
    end do
    end do
    end do
    deallocate(A, B, B_in, B_rf)
end subroutine
subroutine test_ztrsm
    use f77_blas, only: ztrsm, f77_trsm
    use mfi_blas, only: mfi_trsm, mfi_force_gpu, mfi_force_cpu

    integer, parameter :: wp = REAL64
    integer :: M, N
    complex(REAL64), allocatable :: A(:,:), B(:,:), B_in(:,:), B_rf(:,:)
    complex(REAL64) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag
    integer :: ir, ic

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
    call test_get_2d(2000, M)
    N = M

    allocate(A(M,M), B(M,N), B_in(M,N), B_rf(M,N))
    ! Generate a well-conditioned triangular matrix
block
    real(REAL64) :: re(M,M)
    real(REAL64) :: im(M,M)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re(M,N)
    real(REAL64) :: im(M,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL64)
end block
    ! Zero out half, scale off-diagonal, set diagonal to [1,2]
    do ir = 1, M
        do ic = 1, M
            if (ic > ir) then
                A(ir, ic) = 0.0_wp
            else if (ic == ir) then
                A(ir, ic) = 1.0_wp + real(ir, wp) * 0.1_wp
            else
                A(ir, ic) = (A(ir, ic) - 0.5_wp) * 0.2_wp
            end if
        end do
    end do

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side   = sides(i_side)
        uplo   = uplos(i_uplo)
        transa = transas(i_transa)
        diag   = diags(i_diag)

        B_rf = B
        call ztrsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
    call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'mfi_trsm:' // side // uplo // transa // diag // ":&
        & mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
    call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), 'GPU:mfi_trsm:' // side // uplo // transa // diag //&
        & ": mismatch")
#endif

    end do
    end do
    end do
    end do
    deallocate(A, B, B_in, B_rf)
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

