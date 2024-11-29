#:include "common.fpp"
#:include "test/blas/asum_nrm2.fypp"
#:include "test/blas/axpy.fypp"
#:include "test/blas/rot.fypp"
#:include "test/blas/scal.fypp"
program main
use iso_fortran_env
implicit none
$:test_run('?axpy', DEFAULT_TYPES)
$:test_run('?asum', DEFAULT_TYPES, f=MIX_REAL_COMPLEX)
$:test_run('?nrm2', DEFAULT_TYPES, f=MIX_REAL_COMPLEX)
$:test_run('?rot',  DEFAULT_TYPES + COMPLEX_REAL_TYPES)
$:test_run('?scal', DEFAULT_TYPES + COMPLEX_REAL_TYPES)

contains
$:test_implement('?axpy', DEFAULT_TYPES, axpy)
$:test_implement('?asum', DEFAULT_TYPES, asum, f=MIX_REAL_COMPLEX)
$:test_implement('?nrm2', DEFAULT_TYPES, asum, f=MIX_REAL_COMPLEX)
$:test_implement('?rot',  DEFAULT_TYPES, rot)
$:test_implement('?rot',  COMPLEX_TYPES, rot_mixed, MIX_COMPLEX_REAL)
$:test_implement('?scal', DEFAULT_TYPES, scal)
$:test_implement('?scal', COMPLEX_TYPES, scal_mixed, MIX_COMPLEX_REAL)

    pure subroutine assert(test, msg)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        if (.not. test) then
            error stop msg
        end if
    end subroutine
end program
