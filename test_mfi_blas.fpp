#:include "common.fpp"

program test_mfi_blas
    use iso_fortran_env
    use mfi_blas
    use f77_blas
    implicit none
    integer, parameter :: N = 2000
    real(REAL64) :: A(N,N), B(N,N), C(N,N), D(N,N)
    real(REAL64) :: X(N), Y(N), Z(N)
    real(REAL64) :: alpha, beta
    integer :: i, j, k
    C = .0_REAL64
    D = .0_REAL64
    Y = .0_REAL64
    Z = .0_REAL64
    call random_number(A)
    call random_number(B)
    call random_number(X)

    ! BLAS 1
    call test_iamax
    call test_iamin
    ! BLAS 2
    call test_gemv
    ! BLAS 3
    call test_gemm

contains

    subroutine test_iamax
        @:timeit("time f77_iamax: ", { i = f77_iamax(N,X,1) })
        @:timeit("time mfi_iamax: ", { j = mfi_iamax(X)     })
        @:timeit("time maxloc:    ", { k = maxloc(X,1)      })
        call assert(i == j .and. j == k)
    end subroutine

    subroutine test_iamin
        @:timeit("time f77_iamin: ", { i = f77_iamin(N,X,1) })
        @:timeit("time mfi_iamin: ", { j = mfi_iamin(X)     })
        @:timeit("time minloc:    ", { k = minloc(X,1)      })
        call assert(i == j .and. j == k)
    end subroutine

    subroutine test_gemm
        @:timeit("time f77_gemm: ", { call f77_gemm('N', 'N', N, N, N, 1._REAL64, A, N, B, N, 0._REAL64, C, N) })
        @:timeit("time mfi_gemm: ", { call mfi_gemm(A,B,C) })
        @:timeit("time matmul:   ", { D = matmul(A,B)      })
        call assert(all(is_almost_equal(C,D)))

        @:timeit("time f77_gemm:               ", { call f77_gemm('T', 'N', N, N, N, 1._REAL64, A, N, B, N, 0._REAL64, C, N) })
        @:timeit("time mfi_gemm, transa=T:     ", { call mfi_gemm(A,B,C,transa='T') })
        @:timeit("time matmul,   transpose(A): ", { D = matmul(transpose(A),B)      })
        call assert(all(is_almost_equal(C,D)))
    end subroutine

    subroutine test_gemv
        @:timeit("time f77_gemv: ", { call f77_gemv('N', N, N, 1._REAL64, A, N, X, 1, 0._REAL64, Y, 1) })
        @:timeit("time mfi_gemv: ", { call mfi_gemv(A,X,Y) })
        @:timeit("time matmul:   ", { Z = matmul(A,X)      })
        call assert(all(is_almost_equal(Y,Z)))

        @:timeit("time f77_gemv:               ", { call f77_gemv('T', N, N, 1._REAL64, A, N, X, 1, 0._REAL64, Y, 1) })
        @:timeit("time mfi_gemv, transa=T:     ", { call mfi_gemv(A,X,Y,trans='T') })
        @:timeit("time matmul,   transpose(A): ", { Z = matmul(transpose(A),X)     })
        call assert(all(is_almost_equal(Y,Z)))
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
