

program test_axpy
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_saxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against saxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_daxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against daxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_caxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against caxpy", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zaxpy 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_axpy against zaxpy", t2-t1
end block
contains
subroutine test_saxpy
    use f77_blas, only: saxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_saxpy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    real(REAL32) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(REAL32) :: alpha

    call random_number(X)
    call random_number(Y)
    call random_number(alpha)

    x_in = X
    y_in = Y
    call saxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_saxpy(x_in,y_in,alpha)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

end subroutine
subroutine test_daxpy
    use f77_blas, only: daxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_daxpy

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    real(REAL64) :: alpha

    call random_number(X)
    call random_number(Y)
    call random_number(alpha)

    x_in = X
    y_in = Y
    call daxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_daxpy(x_in,y_in,alpha)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

end subroutine
subroutine test_caxpy
    use f77_blas, only: caxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_caxpy

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    complex(REAL32) :: alpha

    real(REAL32) :: rnd_vector(N), rnd
    call random_number(rnd_vector)
    X%re = rnd_vector
    call random_number(rnd_vector)
    X%im = rnd_vector
    call random_number(rnd_vector)
    Y%re = rnd_vector
    call random_number(rnd_vector)
    Y%im = rnd_vector
    call random_number(rnd)
    alpha%re = rnd
    call random_number(rnd)
    alpha%im= rnd

    x_in = X
    y_in = Y
    call caxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_caxpy(x_in,y_in,alpha)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

end subroutine
subroutine test_zaxpy
    use f77_blas, only: zaxpy, f77_axpy
    use mfi_blas, only: mfi_axpy, mfi_zaxpy

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: x(N), Y(N), &
                x_in(N), y_in(N), &
                x_rf(N), y_rf(N)
    complex(REAL64) :: alpha

    real(REAL64) :: rnd_vector(N), rnd
    call random_number(rnd_vector)
    X%re = rnd_vector
    call random_number(rnd_vector)
    X%im = rnd_vector
    call random_number(rnd_vector)
    Y%re = rnd_vector
    call random_number(rnd_vector)
    Y%im = rnd_vector
    call random_number(rnd)
    alpha%re = rnd
    call random_number(rnd)
    alpha%im= rnd

    x_in = X
    y_in = Y
    call zaxpy(N, alpha, x_in, 1, y_in, 1)
    x_rf = x_in
    y_rf = y_in

    x_in = X
    y_in = Y
    call f77_axpy(N, alpha, x_in, 1, y_in, 1)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_zaxpy(x_in,y_in,alpha)
    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

    x_in = X
    y_in = Y
    call mfi_axpy(x_in, y_in, alpha)

    call assert(all(abs(x_in - x_rf) < sqrt(epsilon(1.0_wp))) .and. all(abs(y_in - y_rf) < sqrt(epsilon(1.0_wp))), "different&
        & results")

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

