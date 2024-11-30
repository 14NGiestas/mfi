#:mute
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
#:set COLLECT = [                                             &
    ('?lamch',REAL_TYPES,                       lamch),       &
    ('?dot',  REAL_TYPES,                       dot_product), &
    ('?dotc', COMPLEX_TYPES,                    dot_product), &
    ('?dotu', COMPLEX_TYPES,                    dot_product), &
    ('?copy', DEFAULT_TYPES,                    copy),        &
    ('?swap', DEFAULT_TYPES,                    swap),        &
    ('?axpy', DEFAULT_TYPES,                    axpy),        &
    ('?gemv', DEFAULT_TYPES,                    gemv),        &
    ('?gemm', DEFAULT_TYPES,                    gemm),        &
    ('?asum', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),   &
    ('?nrm2', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),   &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),         &
    ('?rotg', DEFAULT_TYPES, rotg),                           &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),        &
]
#:endmute
program main
use iso_fortran_env
implicit none
#:for name, supported_types, code in COLLECT
$:test_run(name, supported_types)
#:endfor

#:if defined('MFI_EXTENSIONS')
$:test_run('i?amin',DEFAULT_TYPES)
$:test_run('i?amax',DEFAULT_TYPES)
#:endif

contains

#:for name, supported_types, code in COLLECT
$:test_implement(name, supported_types,code)
#:endfor

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
