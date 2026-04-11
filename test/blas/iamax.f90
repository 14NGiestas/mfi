

program test_iamax
use iso_fortran_env
implicit none
#if defined(MFI_EXTENSIONS)
block
real :: t1, t2
call cpu_time(t1)
 call test_isamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against isamax", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_idamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against idamax", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_icamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against icamax", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_izamax 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_iamax against izamax", t2-t1
end block
contains
subroutine test_isamax
    use f77_blas, only: isamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_isamax

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: array(N)
    integer :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = isamax(N, array, 1)
    res(2) = mfi_isamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "sequential array mismatch")

    call random_number(array)
    res(1) = isamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_isamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "random array mismatch")

end subroutine
subroutine test_idamax
    use f77_blas, only: idamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_idamax

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: array(N)
    integer :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = idamax(N, array, 1)
    res(2) = mfi_idamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "sequential array mismatch")

    call random_number(array)
    res(1) = idamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_idamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "random array mismatch")

end subroutine
subroutine test_icamax
    use f77_blas, only: icamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_icamax

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: array(N)
    integer :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = icamax(N, array, 1)
    res(2) = mfi_icamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "sequential array mismatch")

block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im, kind=REAL32)
end block
    res(1) = icamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_icamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "random array mismatch")

end subroutine
subroutine test_izamax
    use f77_blas, only: izamax, f77_iamax
    use mfi_blas, only: mfi_iamax, mfi_izamax

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: array(N)
    integer :: res(4)
    integer :: ii

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

    do ii = 1, N
        array(ii) = real(ii, wp)
    end do
    res(1) = izamax(N, array, 1)
    res(2) = mfi_izamax(array)
    res(3) = f77_iamax(N, array, 1)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "sequential array mismatch")

block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    array = cmplx(re,im, kind=REAL64)
end block
    res(1) = izamax(N, array, 1)
    res(2) = f77_iamax(N, array, 1)
    res(3) = mfi_izamax(array)
    res(4) = mfi_iamax(array)
    call assert(all(res == res(1)), "random array mismatch")

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

#else
    write(*,*) 'i?amax tests skipped: extensions not enabled'
#endif
end program

