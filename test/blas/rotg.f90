

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
    real(REAL32) :: a, b, s
    real(REAL32) :: c

    real(REAL32) :: a_in, b_in, s_in
    real(REAL32) :: c_in

    real(REAL32) :: a_rf, b_rf, s_rf
    real(REAL32) :: c_rf
    integer :: i

    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(s)

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call srotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_drotg
    use f77_blas, only: drotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 200
    real(REAL64) :: a, b, s
    real(REAL64) :: c

    real(REAL64) :: a_in, b_in, s_in
    real(REAL64) :: c_in

    real(REAL64) :: a_rf, b_rf, s_rf
    real(REAL64) :: c_rf
    integer :: i

    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(s)

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call drotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_crotg
    use f77_blas, only: crotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 200
    complex(REAL32) :: a, b, s
    real(REAL32) :: c

    complex(REAL32) :: a_in, b_in, s_in
    real(REAL32) :: c_in

    complex(REAL32) :: a_rf, b_rf, s_rf
    real(REAL32) :: c_rf
    integer :: i

block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    a = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    b = cmplx(re,im)
end block
    call random_number(c)
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    s = cmplx(re,im)
end block

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call crotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
    end do

end subroutine
subroutine test_zrotg
    use f77_blas, only: zrotg
    use mfi_blas, only: mfi_rotg

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 200
    complex(REAL64) :: a, b, s
    real(REAL64) :: c

    complex(REAL64) :: a_in, b_in, s_in
    real(REAL64) :: c_in

    complex(REAL64) :: a_rf, b_rf, s_rf
    real(REAL64) :: c_rf
    integer :: i

block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    a = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    b = cmplx(re,im)
end block
    call random_number(c)
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    s = cmplx(re,im)
end block

    do i=1,N
        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call zrotg(a_in, b_in, c_in, s_in)
        a_rf = a_in
        b_rf = b_in
        c_rf = c_in
        s_rf = s_in

        a_in = a
        b_in = b
        c_in = c
        s_in = s
        call mfi_rotg(a_in, b_in, c_in, s_in)

        call assert(a_in == a_rf .and. &
                    b_in == b_rf .and. &
                    s_in == s_rf .and. &
                    c_in == c_rf, "different results")
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

subroutine report_test_result(test_name, success)
    character(*), intent(in) :: test_name
    logical, intent(in) :: success

    if (success) then
        write(*, '(A, ": ", A)') trim(test_name), 'PASSED'
    else
        write(*, '(A, ": ", A)') trim(test_name), 'FAILED'
    end if
end subroutine

end program

