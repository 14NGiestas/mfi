#:include "common.fpp"

program test_mfi_lapack
    use iso_fortran_env
    use mfi_lapack
    use f77_lapack
    implicit none
    integer, parameter :: N = 2000
    real(REAL64) :: A(N,N)
    real(REAL64) :: X(N), B(N)

contains

    subroutine test_axpy
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
