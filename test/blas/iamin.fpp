
#:include "common.fpp"
#:include "test/blas/macros/iamin_iamax.fypp"

program test_iamin
use iso_fortran_env
implicit none
#:if defined('MFI_EXTENSIONS')
$:test_run('i?amin', DEFAULT_TYPES)
contains
$:test_implement('i?amin', DEFAULT_TYPES, iamin_iamax)

#:include "test/assert.inc"

end program
#:else
! Extensions not enabled, provide minimal program to avoid compilation errors
write(*,*) 'i?amin tests skipped: extensions not enabled'
end program
#:endif

