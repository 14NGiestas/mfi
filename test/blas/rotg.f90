

program test_rotg
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_srotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against srotg", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_drotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against drotg", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_crotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against crotg", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zrotg 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_rotg against zrotg", t2-t1
end block
contains
subroutine test_srotg
    use f77_blas, only: srotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 200
    real(REAL32) :: a, b, s, a_in, b_in, s_in, a_rf, b_rf, s_rf
    real(REAL32) :: c, c_in, c_rf
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
    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(s)

    do i=1,N
        a_in = a; b_in = b; c_in = c; s_in = s
        call srotg(a_in, b_in, c_in, s_in)
        a_rf = a_in; b_rf = b_in; c_rf = c_in; s_rf = s_in

        a_in = a; b_in = b; c_in = c; s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(abs(a_in - a_rf) < sqrt(epsilon(1.0_wp)), "rotg:a mismatch")
        call assert(abs(b_in - b_rf) < sqrt(epsilon(1.0_wp)), "rotg:b mismatch")
        call assert(abs(c_in - c_rf) < sqrt(epsilon(1.0_wp)), "rotg:c mismatch")
        call assert(abs(s_in - s_rf) < sqrt(epsilon(1.0_wp)), "rotg:s mismatch")
    end do

end subroutine
subroutine test_drotg
    use f77_blas, only: drotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 200
    real(REAL64) :: a, b, s, a_in, b_in, s_in, a_rf, b_rf, s_rf
    real(REAL64) :: c, c_in, c_rf
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
    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(s)

    do i=1,N
        a_in = a; b_in = b; c_in = c; s_in = s
        call drotg(a_in, b_in, c_in, s_in)
        a_rf = a_in; b_rf = b_in; c_rf = c_in; s_rf = s_in

        a_in = a; b_in = b; c_in = c; s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(abs(a_in - a_rf) < sqrt(epsilon(1.0_wp)), "rotg:a mismatch")
        call assert(abs(b_in - b_rf) < sqrt(epsilon(1.0_wp)), "rotg:b mismatch")
        call assert(abs(c_in - c_rf) < sqrt(epsilon(1.0_wp)), "rotg:c mismatch")
        call assert(abs(s_in - s_rf) < sqrt(epsilon(1.0_wp)), "rotg:s mismatch")
    end do

end subroutine
subroutine test_crotg
    use f77_blas, only: crotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 200
    complex(REAL32) :: a, b, s, a_in, b_in, s_in, a_rf, b_rf, s_rf
    real(REAL32) :: c, c_in, c_rf
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
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    a = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    b = cmplx(re,im, kind=REAL32)
end block
    call random_number(c)
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    s = cmplx(re,im, kind=REAL32)
end block

    do i=1,N
        a_in = a; b_in = b; c_in = c; s_in = s
        call crotg(a_in, b_in, c_in, s_in)
        a_rf = a_in; b_rf = b_in; c_rf = c_in; s_rf = s_in

        a_in = a; b_in = b; c_in = c; s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(abs(a_in - a_rf) < sqrt(epsilon(1.0_wp)), "rotg:a mismatch")
        call assert(abs(b_in - b_rf) < sqrt(epsilon(1.0_wp)), "rotg:b mismatch")
        call assert(abs(c_in - c_rf) < sqrt(epsilon(1.0_wp)), "rotg:c mismatch")
        call assert(abs(s_in - s_rf) < sqrt(epsilon(1.0_wp)), "rotg:s mismatch")
    end do

end subroutine
subroutine test_zrotg
    use f77_blas, only: zrotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 200
    complex(REAL64) :: a, b, s, a_in, b_in, s_in, a_rf, b_rf, s_rf
    real(REAL64) :: c, c_in, c_rf
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
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    a = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    b = cmplx(re,im, kind=REAL64)
end block
    call random_number(c)
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    s = cmplx(re,im, kind=REAL64)
end block

    do i=1,N
        a_in = a; b_in = b; c_in = c; s_in = s
        call zrotg(a_in, b_in, c_in, s_in)
        a_rf = a_in; b_rf = b_in; c_rf = c_in; s_rf = s_in

        a_in = a; b_in = b; c_in = c; s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(abs(a_in - a_rf) < sqrt(epsilon(1.0_wp)), "rotg:a mismatch")
        call assert(abs(b_in - b_rf) < sqrt(epsilon(1.0_wp)), "rotg:b mismatch")
        call assert(abs(c_in - c_rf) < sqrt(epsilon(1.0_wp)), "rotg:c mismatch")
        call assert(abs(s_in - s_rf) < sqrt(epsilon(1.0_wp)), "rotg:s mismatch")
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

