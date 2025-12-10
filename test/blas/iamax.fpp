
#:include "common.fpp"
#:include "test/blas/macros/iamin_iamax.fypp"

program test_iamax
use iso_fortran_env
implicit none
#:if defined('MFI_EXTENSIONS')
$:test_run('i?amax', DEFAULT_TYPES)
contains
$:test_implement('i?amax', DEFAULT_TYPES, iamin_iamax)

#:include "test/assert.inc"

end program
#:else
! Extensions not enabled, provide minimal program to avoid compilation errors
write(*,*) 'i?amax tests skipped: extensions not enabled'
end program
#:endif

