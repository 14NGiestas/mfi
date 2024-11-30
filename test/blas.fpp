#:include "common.fpp"
#:include "test/blas/asum_nrm2.fypp"
#:include "test/blas/iamin_iamax.fypp"
#:include "test/blas/dot_product.fypp"
#:include "test/blas/axpy.fypp"
#:include "test/blas/rot.fypp"
#:include "test/blas/rotg.fypp"
#:include "test/blas/copy.fypp"
#:include "test/blas/swap.fypp"
#:include "test/blas/scal.fypp"
#:include "test/blas/lamch.fypp"
#:include "test/blas/gemv.fypp"
#:include "test/blas/gemm.fypp"
program main
use iso_fortran_env
implicit none
$:test_run('?lamch',REAL_TYPES)
$:test_run('?dot',  REAL_TYPES)
$:test_run('?dotu', COMPLEX_TYPES)
$:test_run('?dotc', COMPLEX_TYPES)
$:test_run('?copy', DEFAULT_TYPES)
$:test_run('?swap', DEFAULT_TYPES)
$:test_run('?axpy', DEFAULT_TYPES)
$:test_run('?gemv', DEFAULT_TYPES)
$:test_run('?gemm', DEFAULT_TYPES)
$:test_run('?asum', DEFAULT_TYPES, MIX_REAL_COMPLEX)
$:test_run('?nrm2', DEFAULT_TYPES, MIX_REAL_COMPLEX)
$:test_run('?rot',  DEFAULT_TYPES + COMPLEX_REAL_TYPES)
$:test_run('?rotg', DEFAULT_TYPES)
$:test_run('?scal', DEFAULT_TYPES + COMPLEX_REAL_TYPES)

#:if defined('MFI_EXTENSIONS')
$:test_run('i?amin',DEFAULT_TYPES)
$:test_run('i?amax',DEFAULT_TYPES)
#:endif

contains

$:test_implement('?lamch',REAL_TYPES, lamch)
$:test_implement('?dot',  REAL_TYPES,    dot_product)
$:test_implement('?dotc', COMPLEX_TYPES, dot_product)
$:test_implement('?dotu', COMPLEX_TYPES, dot_product)
$:test_implement('?copy', DEFAULT_TYPES, copy)
$:test_implement('?swap', DEFAULT_TYPES, swap)
$:test_implement('?axpy', DEFAULT_TYPES, axpy)
$:test_implement('?gemv', DEFAULT_TYPES, gemv)
$:test_implement('?gemm', DEFAULT_TYPES, gemm)
$:test_implement('?asum', DEFAULT_TYPES, asum_nrm2, MIX_REAL_COMPLEX)
$:test_implement('?nrm2', DEFAULT_TYPES, asum_nrm2, MIX_REAL_COMPLEX)
$:test_implement('?rot',  DEFAULT_TYPES, rot)
$:test_implement('?rot',  COMPLEX_TYPES, rot_mixed, MIX_COMPLEX_REAL)
$:test_implement('?rotg', DEFAULT_TYPES, rotg)
$:test_implement('?scal', DEFAULT_TYPES, scal)
$:test_implement('?scal', COMPLEX_TYPES, scal_mixed, MIX_COMPLEX_REAL)

#:if defined('MFI_EXTENSIONS')
$:test_implement('i?amin',DEFAULT_TYPES, iamin_iamax)
$:test_implement('i?amax',DEFAULT_TYPES, iamin_iamax)
#:endif

    pure subroutine assert(test, msg)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        if (.not. test) then
            error stop msg
        end if
    end subroutine
end program
