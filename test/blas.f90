
program test_mfi_blas
    use iso_fortran_env
    use mfi_blas
    use f77_blas
    implicit none
    integer, parameter :: N = 2000
    real(REAL64) :: A(N,N), B(N,N), C(N,N), D(N,N)
    real(REAL64) :: X(N), Y(N), Z(N), U(N)
    real(REAL64) :: alpha, beta
    integer :: i, j, k

    ! Auxiliary
    call test_lamch
    ! BLAS 1
    call test_axpy
    call test_copy
    call test_rotmg
    call test_swap
    call test_iamax
    call test_iamin
    ! BLAS 2
    call test_gemv
    ! BLAS 3
    call test_gemm

contains
    subroutine test_defaults
        C = .0_REAL64
        D = .0_REAL64
        Y = .0_REAL64
        Z = .0_REAL64
        alpha = 1.0_REAL64
        beta  = 0.0_REAL64
        call random_number(A)
        call random_number(B)
        call random_number(X)
    end subroutine

    subroutine test_rotmg
        real(REAL64) :: x1, y1, d1, d2, params(5)
        real(REAL64) :: expected(5)
        expected = [-1.d0, 1.6110934624105326d-6, -2.44140625d-4, 2.44140625d-4, 1.62760416d-6]
        d1 = 5.9d-8; d2 = 5.960464d-8; x1 = 1.d0; y1 = 150.d0; params = .0_REAL64
block
real :: t1, t2
call cpu_time(t1)
 call f77_rotmg(d1, d2, x1, y1, params) 
call cpu_time(t2)
print '(A,G0)', "time f77_rotmg: ", t2-t1
end block
        call assert(all(is_almost_equal(params, expected)))
        d1 = 5.9d-8; d2 = 5.960464d-8; x1 = 1.d0; y1 = 150.d0; params = .0_REAL64
block
real :: t1, t2
call cpu_time(t1)
 call mfi_rotmg(d1, d2, x1, y1, params) 
call cpu_time(t2)
print '(A,G0)', "time mfi_rotmg: ", t2-t1
end block
        call assert(all(is_almost_equal(params, expected)))
        d1 = 5.9d-8; d2 = 5.960464d-8; x1 = 1.d0; y1 = 150.d0; params = .0_REAL64
block
real :: t1, t2
call cpu_time(t1)
 call drotmg(d1, d2, x1, y1, params) 
call cpu_time(t2)
print '(A,G0)', "time drotmg:    ", t2-t1
end block
        call assert(all(is_almost_equal(params, expected)))
    end subroutine

    subroutine test_axpy
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 call f77_axpy(N,2*alpha,X,1,Y,1) 
call cpu_time(t2)
print '(A,G0)', "time f77_axpy: ", t2-t1
end block
        Y = .0_REAL64
block
real :: t1, t2
call cpu_time(t1)
 call mfi_axpy(X,Y,2*alpha)       
call cpu_time(t2)
print '(A,G0)', "time mfi_axpy: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 Z=2*X                            
call cpu_time(t2)
print '(A,G0)', "time fortran:  ", t2-t1
end block
        call assert(all(is_almost_equal(2*X,Y) .and. is_almost_equal(2*X,Z)))
    end subroutine

    subroutine test_copy
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 call f77_copy(N,X,1,Y,1) 
call cpu_time(t2)
print '(A,G0)', "time f77_copy: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call mfi_copy(X,Y)       
call cpu_time(t2)
print '(A,G0)', "time mfi_copy: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 Z=X                      
call cpu_time(t2)
print '(A,G0)', "time fortran:  ", t2-t1
end block
        call assert(all(is_almost_equal(X,Y) .and. is_almost_equal(X,Z)))
    end subroutine

    subroutine test_swap
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 call f77_swap(N,X,1,Y,1) 
call cpu_time(t2)
print '(A,G0)', "time f77_swap: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call mfi_swap(X,Y)       
call cpu_time(t2)
print '(A,G0)', "time mfi_swap: ", t2-t1
end block
    end subroutine

    subroutine test_iamax
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 i = f77_iamax(N,X,1) 
call cpu_time(t2)
print '(A,G0)', "time f77_iamax: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 j = mfi_iamax(X)     
call cpu_time(t2)
print '(A,G0)', "time mfi_iamax: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 k = maxloc(X,1)      
call cpu_time(t2)
print '(A,G0)', "time maxloc:    ", t2-t1
end block
        call assert(i == j .and. j == k)
    end subroutine

    subroutine test_iamin
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 i = f77_iamin(N,X,1) 
call cpu_time(t2)
print '(A,G0)', "time f77_iamin: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 j = mfi_iamin(X)     
call cpu_time(t2)
print '(A,G0)', "time mfi_iamin: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 k = minloc(X,1)      
call cpu_time(t2)
print '(A,G0)', "time minloc:    ", t2-t1
end block
        call assert(i == j .and. j == k)
    end subroutine

    subroutine test_gemm
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 call f77_gemm('N', 'N', N, N, N, alpha, A, N, B, N, beta, C, N) 
call cpu_time(t2)
print '(A,G0)', "time f77_gemm: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call mfi_gemm(A,B,C) 
call cpu_time(t2)
print '(A,G0)', "time mfi_gemm: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 D = matmul(A,B)      
call cpu_time(t2)
print '(A,G0)', "time matmul:   ", t2-t1
end block
        call assert(all(is_almost_equal(C,D)))

block
real :: t1, t2
call cpu_time(t1)
 call f77_gemm('T', 'N', N, N, N, alpha, A, N, B, N, beta, C, N) 
call cpu_time(t2)
print '(A,G0)', "time f77_gemm, transa=T:     ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call mfi_gemm(A,B,C,transa='T') 
call cpu_time(t2)
print '(A,G0)', "time mfi_gemm, transa=T:     ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 D = matmul(transpose(A),B)      
call cpu_time(t2)
print '(A,G0)', "time matmul,   transpose(A): ", t2-t1
end block
        call assert(all(is_almost_equal(C,D)))
    end subroutine

    subroutine test_gemv
        call test_defaults
block
real :: t1, t2
call cpu_time(t1)
 call f77_gemv('N', N, N, alpha, A, N, X, 1, beta, Y, 1) 
call cpu_time(t2)
print '(A,G0)', "time f77_gemv: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call mfi_gemv(A,X,Y) 
call cpu_time(t2)
print '(A,G0)', "time mfi_gemv: ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 Z = matmul(A,X)      
call cpu_time(t2)
print '(A,G0)', "time matmul:   ", t2-t1
end block
        call assert(all(is_almost_equal(Y,Z)))

block
real :: t1, t2
call cpu_time(t1)
 call f77_gemv('T', N, N, alpha, A, N, X, 1, beta, Y, 1) 
call cpu_time(t2)
print '(A,G0)', "time f77_gemv: trans=T:      ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call mfi_gemv(A,X,Y,trans='T') 
call cpu_time(t2)
print '(A,G0)', "time mfi_gemv, trans=T:      ", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 Z = matmul(transpose(A),X)     
call cpu_time(t2)
print '(A,G0)', "time matmul,   transpose(A): ", t2-t1
end block
        call assert(all(is_almost_equal(Y,Z)))
    end subroutine

    subroutine test_lamch
        real(REAL32) :: sa
        real(REAL64) :: da
        sa = mfi_lamch('E',1.0_REAL32)
        da = mfi_lamch('E',1.0_REAL64)
        call assert(sa > da)
    end subroutine

    pure subroutine assert(test)
        logical, intent(in) :: test
        if (.not. test) then
            error stop 'different results'
        end if
    end subroutine

    logical pure elemental function is_almost_equal(x, y)
        real(REAL64), intent(in) :: x, y
        is_almost_equal = abs(x-y) < 10**6*epsilon(x)
    end function

end program
