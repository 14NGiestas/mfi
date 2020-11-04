program test_mfi_blas
    use iso_fortran_env
    use mfi_blas
    implicit none
    integer, parameter :: N = 1000
    real(REAL64) :: t1, t2
    real(REAL64) :: A(N,N), B(N,N), C(N,N), D(N,N)
    real(REAL64) :: X(N), Y(N), Z(N)
    real(REAL64) :: alpha, beta
    integer :: i, j
    C = .0_REAL64
    D = .0_REAL64
    Y = .0_REAL64
    Z = .0_REAL64
    call random_number(A)
    call random_number(B)
    call random_number(X)

    call test_iamax
    call test_iamin
    call test_gemm
    call test_gemv

contains

    subroutine test_iamax
        print '("test_iamax")'
        call cpu_time(t1)
        i = iamax(X)
        call cpu_time(t2)
        print '("time iamax: ",G0)', t2-t1

        call cpu_time(t1)
        j = maxloc(X,1)
        call cpu_time(t2)
        print '("time maxloc: ",G0)', t2-t1

        call assert(i == j)
    end subroutine

    subroutine test_iamin
        print '("test_iamin")'
        call cpu_time(t1)
        i = iamin(X)
        call cpu_time(t2)
        print '("time iamin: ",G0)', t2-t1

        call cpu_time(t1)
        j = minloc(X,1)
        call cpu_time(t2)
        print '("time minloc: ",G0)', t2-t1

        call assert(i == j)
    end subroutine

    subroutine test_gemm
        print '("test_gemm")'
        call cpu_time(t1)
        call gemm(A,B,C)
        call cpu_time(t2)
        print '("time gemm: ",G0)', t2-t1

        call cpu_time(t1)
        D = matmul(A,B)
        call cpu_time(t2)
        print '("time matmul: ",G0)', t2-t1

        call assert(all(abs(C-D) < 10**6*epsilon(C-D)))
    end subroutine

    subroutine test_gemv
        print '("test_gemv")'
        call cpu_time(t1)
        call gemv(A,X,Y)
        call cpu_time(t2)
        print '("time gemv: ",G0)', t2-t1

        call cpu_time(t1)
        Z = matmul(A,X)
        call cpu_time(t2)
        print '("time matmul: ",G0)', t2-t1

        call assert(all(abs(Z-Y) < 10**6*epsilon(Z-Y)))
    end subroutine

    pure subroutine assert(test)
        logical, intent(in) :: test
        if (.not. test) then
            error stop 'different results'
        end if
    end subroutine

end program
