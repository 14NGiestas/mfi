

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
    real(REAL32) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    real(REAL32) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

    call random_number(M)
    call random_number(X)
    call random_number(Y)
    call random_number(alpha)
    call random_number(beta)


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call sgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_sgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_dgemv
    use f77_blas, only: dgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_dgemv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    real(REAL64) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    real(REAL64) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

    call random_number(M)
    call random_number(X)
    call random_number(Y)
    call random_number(alpha)
    call random_number(beta)


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call dgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_dgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_cgemv
    use f77_blas, only: cgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_cgemv

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    complex(REAL32) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    complex(REAL32) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL32) :: re(N)
    real(REAL32) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im)
end block


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call cgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_cgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
    end do

end subroutine
subroutine test_zgemv
    use f77_blas, only: zgemv, f77_gemv
    use mfi_blas, only: mfi_gemv, mfi_zgemv

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    complex(REAL64) :: M(N,N),    X(N),    Y(N),   &
                M_in(N,N), X_in(N), Y_in(N),&
                M_rf(N,N), X_rf(N), Y_rf(N)
    complex(REAL64) :: alpha, beta
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    character :: trans
    integer :: i

block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    M = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    X = cmplx(re,im)
end block
block
    real(REAL64) :: re(N)
    real(REAL64) :: im(N)
    call random_number(im)
    call random_number(re)
    Y = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    beta = cmplx(re,im)
end block


    do i=1,size(options)
        trans = options(i)

        M_in = M
        X_in = X
        Y_in = Y
        call zgemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        M_rf = M_in
        X_rf = X_in
        Y_rf = Y_in

        M_in = M
        X_in = X
        Y_in = Y
        call f77_gemv(trans, N, N, alpha, M_in, N, X_in, 1, beta, Y_in, 1)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_zgemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)
        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")

        M_in = M
        X_in = X
        Y_in = Y
        call mfi_gemv(M_in,X_in,Y_in,alpha=alpha, beta=beta, trans=trans)

        call assert(all(M_in == M_rf) .and. &
                    all(X_in == X_rf) .and. &
                    all(Y_in == Y_rf), "different results")
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

