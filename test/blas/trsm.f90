

program test_trsm
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_strsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm (CPU) against strsm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dtrsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm (CPU) against dtrsm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ctrsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm (CPU) against ctrsm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ztrsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm (CPU) against ztrsm", t2-t1
end block
contains
subroutine test_strsm
    use f77_blas, only: strsm, f77_trsm
    use mfi_blas, only: mfi_trsm

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 10, N = 10
    real(REAL32) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    real(REAL32) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

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

    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do

    call random_number(B)
    call random_number(alpha)

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
call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL32)), "${mfi}$:${side}${uplo}${transa}${diag}: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:${side}${uplo}${transa}${diag}: mismatch")
#endif

    end do
    end do
    end do
    end do
end subroutine
subroutine test_dtrsm
    use f77_blas, only: dtrsm, f77_trsm
    use mfi_blas, only: mfi_trsm

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 10, N = 10
    real(REAL64) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    real(REAL64) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

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

    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do

    call random_number(B)
    call random_number(alpha)

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
call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL64)), "${mfi}$:${side}${uplo}${transa}${diag}: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
call assert(maxval(abs(B_in - B_rf)) < sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:${side}${uplo}${transa}${diag}: mismatch")
#endif

    end do
    end do
    end do
    end do
end subroutine
subroutine test_ctrsm
    use f77_blas, only: ctrsm, f77_trsm
    use mfi_blas, only: mfi_trsm

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 10, N = 10
    complex(REAL32) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    complex(REAL32) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

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

    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do

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
call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "${mfi}$:${side}${uplo}${transa}${diag}: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL32)), "GPU:${mfi}$:${side}${uplo}${transa}${diag}: mismatch")
#endif

    end do
    end do
    end do
    end do
end subroutine
subroutine test_ztrsm
    use f77_blas, only: ztrsm, f77_trsm
    use mfi_blas, only: mfi_trsm

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 10, N = 10
    complex(REAL64) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    complex(REAL64) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

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

    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do

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
call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "${mfi}$:${side}${uplo}${transa}${diag}: mismatch")

#if defined(MFI_CUBLAS)
        call mfi_force_gpu()
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        call mfi_force_cpu()
call assert(maxval(abs(B_in - B_rf)) < 2.0 * sqrt(epsilon(1.0_REAL64)), "GPU:${mfi}$:${side}${uplo}${transa}${diag}: mismatch")
#endif

    end do
    end do
    end do
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
