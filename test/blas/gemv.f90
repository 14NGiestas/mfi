

program test_gemv
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against sgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against dgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against cgemv", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zgemv 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemv against zgemv", t2-t1
end block
contains
subroutine test_sgemv
    use f77_blas, only: sgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_sgemv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: M(N,N), x(N), y(N), &
                M_in(N,N), x_in(N), y_in(N), &
                M_rf(N,N), x_rf(N), y_rf(N)
    real(REAL32) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

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
    call random_number(M)
    call random_number(x)
    call random_number(y)
    call random_number(alpha)
    call random_number(beta)

    do i=1,size(options)
        trans = options(i)

        M_in = M; x_in = x; y_in = y
        call sgemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
        M_rf = M_in; x_rf = x_in; y_rf = y_in

        M_in = M; x_in = x; y_in = y
        call f77_gemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_sgemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
        call mfi_force_cpu()
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:y: mismatch")
#endif
    end do

end subroutine
subroutine test_dgemv
    use f77_blas, only: dgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_dgemv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: M(N,N), x(N), y(N), &
                M_in(N,N), x_in(N), y_in(N), &
                M_rf(N,N), x_rf(N), y_rf(N)
    real(REAL64) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

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
    call random_number(M)
    call random_number(x)
    call random_number(y)
    call random_number(alpha)
    call random_number(beta)

    do i=1,size(options)
        trans = options(i)

        M_in = M; x_in = x; y_in = y
        call dgemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
        M_rf = M_in; x_rf = x_in; y_rf = y_in

        M_in = M; x_in = x; y_in = y
        call f77_gemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_dgemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:y: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
        call mfi_force_cpu()
call assert(maxval(abs(M_in - M_rf)) < sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:y: mismatch")
#endif
    end do

end subroutine
subroutine test_cgemv
    use f77_blas, only: cgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_cgemv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: M(N,N), x(N), y(N), &
                M_in(N,N), x_in(N), y_in(N), &
                M_rf(N,N), x_rf(N), y_rf(N)
    complex(REAL32) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

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
block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    x = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im, kind=REAL32)
end block

    do i=1,size(options)
        trans = options(i)

        M_in = M; x_in = x; y_in = y
        call cgemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
        M_rf = M_in; x_rf = x_in; y_rf = y_in

        M_in = M; x_in = x; y_in = y
        call f77_gemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${f90}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_cgemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "mfi_${f77}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:y: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
        call mfi_force_cpu()
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:y: mismatch")
#endif
    end do

end subroutine
subroutine test_zgemv
    use f77_blas, only: zgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_zgemv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: M(N,N), x(N), y(N), &
                M_in(N,N), x_in(N), y_in(N), &
                M_rf(N,N), x_rf(N), y_rf(N)
    complex(REAL64) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

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
block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    x = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    y = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im, kind=REAL64)
end block

    do i=1,size(options)
        trans = options(i)

        M_in = M; x_in = x; y_in = y
        call zgemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
        M_rf = M_in; x_rf = x_in; y_rf = y_in

        M_in = M; x_in = x; y_in = y
        call f77_gemv(trans, N, N, alpha, M_in, N, x_in, 1, beta, y_in, 1)
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${f90}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_zgemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "mfi_${f77}$:y: mismatch")

        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:y: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        M_in = M; x_in = x; y_in = y
        call mfi_gemv(M_in, x_in, y_in, alpha=alpha, beta=beta, trans=trans)
        call mfi_force_cpu()
call assert(maxval(abs(M_in - M_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:M: mismatch")
call assert(maxval(abs(x_in - x_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:x: mismatch")
call assert(maxval(abs(y_in - y_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:y: mismatch")
#endif
    end do

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

