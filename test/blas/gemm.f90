

program test_gemm
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_sgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against sgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against dgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against cgemm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zgemm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_gemm against zgemm", t2-t1
end block
contains
subroutine test_sgemm
    use f77_blas, only: sgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_sgemm

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    real(REAL32) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    real(REAL32) :: alpha, beta
    character :: transa, transb
    integer :: i, j

    call random_number(A)
    call random_number(B)
    call random_number(C)
    call random_number(alpha)
    call random_number(beta)

    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call sgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_sgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_dgemm
    use f77_blas, only: dgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_dgemm

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    real(REAL64) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    real(REAL64) :: alpha, beta
    character :: transa, transb
    integer :: i, j

    call random_number(A)
    call random_number(B)
    call random_number(C)
    call random_number(alpha)
    call random_number(beta)

    do i=1,size(options)
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call dgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_dgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_cgemm
    use f77_blas, only: cgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_cgemm

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    complex(REAL32) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    complex(REAL32) :: alpha, beta
    character :: transa, transb
    integer :: i, j

block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im)
end block
block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im)
end block
block
    real(REAL32) :: re(N,N)
    real(REAL32) :: im(N,N)
    call random_number(im)
    call random_number(re)
    C = cmplx(re,im)
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
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call cgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_cgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
    end do
    end do

end subroutine
subroutine test_zgemm
    use f77_blas, only: zgemm, f77_gemm
    use mfi_blas, only: mfi_gemm, mfi_zgemm

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 20
    character, parameter :: options(*) = ['N','n','T','t','C','c']
    complex(REAL64) :: A(N,N),    B(N,N),    C(N,N),   &
                A_in(N,N), B_in(N,N), C_in(N,N),&
                A_rf(N,N), B_rf(N,N), C_rf(N,N)
    complex(REAL64) :: alpha, beta
    character :: transa, transb
    integer :: i, j

block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    A = cmplx(re,im)
end block
block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im)
end block
block
    real(REAL64) :: re(N,N)
    real(REAL64) :: im(N,N)
    call random_number(im)
    call random_number(re)
    C = cmplx(re,im)
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
    do j=1,size(options)
        transa = options(i)
        transb = options(j)

        A_in = A
        B_in = B
        C_in = C
        call zgemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        A_rf = A_in
        B_rf = B_in
        C_rf = C_in

        A_in = A
        B_in = B
        C_in = C
        call f77_gemm(transa, transb, N, N, N, alpha, A_in, N, B_in, N, beta, C_in, N)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_zgemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)
        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")

        A_in = A
        B_in = B
        C_in = C
        call mfi_gemm(A_in,B_in,C_in,alpha=alpha, beta=beta, transa=transa, transb=transb)

        call assert(all(A_in == A_rf) .and. &
                    all(B_in == B_rf) .and. &
                    all(C_in == C_rf), "different results")
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

