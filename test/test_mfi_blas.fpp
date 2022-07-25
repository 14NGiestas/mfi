#:include "common.fpp"

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

    ! BLAS 1
    call test_axpy
    call test_copy
    call test_rotmg
    call test_swap
    call test_iamax
#:if defined('MFI_EXTENSIONS')
    call test_iamin
#:endif
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
        @:timeit("time f77_rotmg: ", { call f77_rotmg(d1, d2, x1, y1, params) })
        call assert(all(is_almost_equal(params, expected)))
        d1 = 5.9d-8; d2 = 5.960464d-8; x1 = 1.d0; y1 = 150.d0; params = .0_REAL64
        @:timeit("time mfi_rotmg: ", { call mfi_rotmg(d1, d2, x1, y1, params) })
        call assert(all(is_almost_equal(params, expected)))
        d1 = 5.9d-8; d2 = 5.960464d-8; x1 = 1.d0; y1 = 150.d0; params = .0_REAL64
        @:timeit("time drotmg:    ", { call drotmg(d1, d2, x1, y1, params) })
        call assert(all(is_almost_equal(params, expected)))
    end subroutine

    subroutine test_axpy
        call test_defaults
        @:timeit("time f77_axpy: ", { call f77_axpy(N,2*alpha,X,1,Y,1) })
        Y = .0_REAL64
        @:timeit("time mfi_axpy: ", { call mfi_axpy(X,Y,2*alpha)       })
        @:timeit("time fortran:  ", { Z=2*X                            })
        call assert(all(is_almost_equal(2*X,Y) .and. is_almost_equal(2*X,Z)))
    end subroutine

    subroutine test_copy
        call test_defaults
        @:timeit("time f77_copy: ", { call f77_copy(N,X,1,Y,1) })
        @:timeit("time mfi_copy: ", { call mfi_copy(X,Y)       })
        @:timeit("time fortran:  ", { Z=X                      })
        call assert(all(is_almost_equal(X,Y) .and. is_almost_equal(X,Z)))
    end subroutine

    subroutine test_swap
        call test_defaults
        @:timeit("time f77_swap: ", { call f77_swap(N,X,1,Y,1) })
        @:timeit("time mfi_swap: ", { call mfi_swap(X,Y)       })
    end subroutine

    subroutine test_iamax
        call test_defaults
        @:timeit("time f77_iamax: ", { i = f77_iamax(N,X,1) })
        @:timeit("time mfi_iamax: ", { j = mfi_iamax(X)     })
        @:timeit("time maxloc:    ", { k = maxloc(X,1)      })
        call assert(i == j .and. j == k)
    end subroutine

#:if defined('MFI_EXTENSIONS')
    subroutine test_iamin
        call test_defaults
        @:timeit("time f77_iamin: ", { i = f77_iamin(N,X,1) })
        @:timeit("time mfi_iamin: ", { j = mfi_iamin(X)     })
        @:timeit("time minloc:    ", { k = minloc(X,1)      })
        call assert(i == j .and. j == k)
    end subroutine
#:endif

    subroutine test_gemm
        call test_defaults
        @:timeit("time f77_gemm: ", { call f77_gemm('N', 'N', N, N, N, alpha, A, N, B, N, beta, C, N) })
        @:timeit("time mfi_gemm: ", { call mfi_gemm(A,B,C) })
        @:timeit("time matmul:   ", { D = matmul(A,B)      })
        call assert(all(is_almost_equal(C,D)))

        @:timeit("time f77_gemm, transa=T:     ", { call f77_gemm('T', 'N', N, N, N, alpha, A, N, B, N, beta, C, N) })
        @:timeit("time mfi_gemm, transa=T:     ", { call mfi_gemm(A,B,C,transa='T') })
        @:timeit("time matmul,   transpose(A): ", { D = matmul(transpose(A),B)      })
        call assert(all(is_almost_equal(C,D)))
    end subroutine

    subroutine test_gemv
        call test_defaults
        @:timeit("time f77_gemv: ", { call f77_gemv('N', N, N, alpha, A, N, X, 1, beta, Y, 1) })
        @:timeit("time mfi_gemv: ", { call mfi_gemv(A,X,Y) })
        @:timeit("time matmul:   ", { Z = matmul(A,X)      })
        call assert(all(is_almost_equal(Y,Z)))

        @:timeit("time f77_gemv: trans=T:      ", { call f77_gemv('T', N, N, alpha, A, N, X, 1, beta, Y, 1) })
        @:timeit("time mfi_gemv, trans=T:      ", { call mfi_gemv(A,X,Y,trans='T') })
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
