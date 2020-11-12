#:include "common.fpp"

program test_mfi_lapack
    use iso_fortran_env
    use mfi_lapack
    use f77_lapack
    implicit none
    integer, parameter :: N = 2000
    real(REAL64) :: S(N,N)
    integer :: i, j, info

    call test_potrf
    call test_potri

contains

    subroutine symmetric
        call random_number(S)
        do j=1,N
            do i=j+1,N
                S(i,j) = S(j,i)
            end do
        end do
    end subroutine

    subroutine positive_definite
        call symmetric
        S = matmul(S,S)
    end subroutine

    subroutine test_potrf
        call positive_definite
        @:timeit("f77_potrf: ", { call f77_potrf('U',N,S,N,info) })
        call assert(info == 0)

        call positive_definite
        @:timeit("mfi_potrf: ", { call mfi_potrf(S,info)         })
        call assert(info == 0)
    end subroutine

    subroutine test_potri
        call positive_definite
        call f77_potrf('U',N,S,N,info)
        @:timeit("f77_potri: ", { call f77_potri('U',N,S,N,info) })
        call assert(info == 0)

        call positive_definite
        call mfi_potrf(S,info)
        @:timeit("mfi_potri: ", { call mfi_potri(S,info)         })
        call assert(info == 0)
    end subroutine

    pure subroutine assert(test)
        logical, intent(in) :: test
        if (.not. test) then
            error stop 'assertion failed'
        end if
    end subroutine

    logical pure elemental function is_almost_equal(x, y)
        real(REAL64), intent(in) :: x, y
        is_almost_equal = abs(x-y) < 10**6*epsilon(x)
    end function

end program
