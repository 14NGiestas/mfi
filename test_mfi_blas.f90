program test_mfi_blas
    use iso_fortran_env
    use mfi_blas
    implicit none
    integer, parameter :: N = 1000
    real(REAL64) :: t1, t2
    real(REAL64) :: A(N,N), B(N,N), C(N,N), D(N,N)
    real(REAL64) :: X(N), Y(N), Z(N)
    call random_number(A)
    call random_number(B)
    call random_number(X)

    call test_gemm
    call test_gemv

contains

    subroutine test_gemm
        call cpu_time(t1)
        call gemm(A,B,C)
        call cpu_time(t2)
        print '("time gemm: ",G0)', t2-t1

        call cpu_time(t1)
        D = matmul(A,B)
        call cpu_time(t2)
        print '("time matmul: ",G0)', t2-t1

        if (all(abs(C-D) < 10**6*epsilon(C-D))) then
            print '("both produce the same results")'
        else
            error stop 'different results'
        end if
    end subroutine

    subroutine test_gemv
        call cpu_time(t1)
        call gemv(A,X,Y)
        call cpu_time(t2)
        print '("time gemv: ",G0)', t2-t1

        call cpu_time(t1)
        Z = matmul(A,X)
        call cpu_time(t2)
        print '("time matmul: ",G0)', t2-t1

        if (all(abs(Z-Y) < 10**6*epsilon(Z-Y))) then
            print '("both produce the same results")'
        else
            error stop 'different results'
        end if
    end subroutine

end program
