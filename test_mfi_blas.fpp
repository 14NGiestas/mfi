#:include "common.fpp"

program test_mfi_blas
    use iso_fortran_env
    use mfi_blas
    implicit none
    integer, parameter :: N = 1000
    real(REAL64) :: t1, t2
    real(REAL64) :: A(N,N), B(N,N), C(N,N), D(N,N)
    real(REAL64) :: X(N), Y(N), Z(N)
    integer :: i, j, k
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
        integer, external :: idamax
        @:timeit("time izmax:  ", { i = idamax(N,X,1) })
        @:timeit("time imax:   ", { j = iamax(X)      })
        @:timeit("time maxloc: ", { k = maxloc(X,1)   })
        call assert(i == j .and. j == k)
    end subroutine

    subroutine test_iamin
        integer, external :: idamin
        @:timeit("time izamin: ", { i = idamin(N,X,1) })
        @:timeit("time iamin:  ", { j = iamin(X)      })
        @:timeit("time minloc: ", { k = minloc(X,1)   })
        call assert(i == j .and. j == k)
    end subroutine

    subroutine test_gemm
        @:timeit("time gemm:   ", { call gemm(A,B,C) })
        @:timeit("time matmul: ", { D = matmul(A,B)  })
        call assert(all(is_almost_equal(C,D)))

        @:timeit("time gemm,   transa=T:     ", { call gemm(A,B,C,transa='T') })
        @:timeit("time matmul, transpose(A): ", { D = matmul(transpose(A),B)  })
        call assert(all(is_almost_equal(C,D)))
    end subroutine

    subroutine test_gemv
        @:timeit("time gemv:   ", { call gemv(A,X,Y) })
        @:timeit("time matmul: ", { Z = matmul(A,X)  })
        call assert(all(is_almost_equal(Y,Z)))

        @:timeit("time gemv,   transa=T:     ", { call gemv(A,X,Y,trans='T') })
        @:timeit("time matmul, transpose(A): ", { Z = matmul(transpose(A),X)  })
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
