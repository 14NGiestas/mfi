
#:include "common.fpp"

program test_iamax
use iso_fortran_env
implicit none
#:if defined('MFI_EXTENSIONS')
#:include "test/blas/macros/iamin_iamax.fypp"
$:test_run('i?amax', DEFAULT_TYPES)
contains
$:test_implement('i?amax', DEFAULT_TYPES, iamin_iamax)

#:include "test/assert.inc"

end program
#:else
    write(*,*) 'i?amax tests skipped: extensions not enabled'
end program
#:endif

